#ifndef __PRECONDITIONER
#define __PRECONDITIONER

// Written by: C. Nguyen & P. Fernandez

static Int indexofSmallestElement(double* array, Int size)
{
    Int index = 0;
    
    for(Int i = 1; i < size; i++)
    {
        if (isnan(array[i]) || isinf(array[i]))
            error("NaN or Inf detected in indexofSmallestElement.\n");
        
        if(array[i] < array[index])
            index = i;
    }
    return index;
}

static void computeNoOrdering(sysstruct &sys, Int preconditioner)
{
    Int nrows;
    
    if (preconditioner == 0)
        nrows = sys.nrowsP;
    else if (preconditioner == 1)
        nrows = sys.nrowsP;
    else if (preconditioner == -1 || preconditioner == 2)
        nrows = 0;
    else
        error("Preconditioner flag not recognized in computeNoOrdering.\n");
    
    for (Int i=0; i<nrows; i++) {
        sys.unordered2ordered[i] = i;
        sys.ordered2unordered[i] = i;
    }
}

void computeApproximateMDFordering(sysstruct &sys, Int preconditioner)
{
    /* TODO: Further optimization: If necessary, this function can be optimized by not recomputing w[i] when needs to be updated (we can just remove the contributions that no longer exist). */

    /* Note:
     * In this function, entity refers to faces in HDG and to face nodes in EDG and IEDG.
     * We assume all blocks are of the same size (i.e., nch (EDG) or nch*npf (UDH) are constant throughout the mesh).
     *      If not, an array indicating the size of each block and an array indicating the starting position each block in the matrix J would be required (easy fix).
     * The complexity of the MDF ordering is O(N^2) due to the indexofSmallestElement operation.
     * */

    char chn = 'N';
    double one = 1.0, zero = 0.0;
    
    Int nrows, nblks;
    Int *rowpts, *colind;
    if (preconditioner == 0) {
        nrows = sys.numEntities;
        nblks = sys.numBlocks;
        rowpts = &sys.ent2entStart[0];
        colind = &sys.ent2ent[0];
    }
    else if (preconditioner == 1) {
        nrows = sys.nrowsP;//BJ_nrows;
        nblks = sys.numBlocks;
        rowpts = &sys.ent2entStart[0];
        colind = &sys.ent2ent[0];
        //nrows = sys.BJ_nrows;
        //nblks = sys.BJ_nblks;
        //rowpts = &sys.BJ_rowpts[0];
        //colind = &sys.BJ_colind[0];
    }
    else
        error("Preconditioner flag not recognized in computeApproximateMDFordering.\n");
    
    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int inc = 1, i,j,k,m,n, info;
    Int rk, kk, nrk, ri, rj, nri, pij, pki, ii, ik, kj, rl, pl, pim, nrl, plk, shared;
    double Cik, Dik, Ckj, bign = 1.0e100;
    
    vector<vector<Int> >* ipivs = &sys.ipivs;
    Int lwork = bsz;
    vector<vector<double> >* works = &sys.works;
    Int *pl_local = new Int[sys.noThreads];
    
//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for private(info, j, ii, nri, pij)
        for (Int ri=0; ri<nrows; ri++) {
            ii = rowpts[ri];
            nri = rowpts[ri+1] - rowpts[ri];

            /* HiiInv <- inv(Hg_ii) */
            DCOPY(&bsz2, &sys.Hg[bsz2*ii], &inc, &sys.HiiInvs[this_thread][0], &inc);
            DGETRF(&bsz, &bsz, &sys.HiiInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &info);
            DGETRI(&bsz, &sys.HiiInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
            
            for (j=0; j<nri; j++) {
                pij = ii+j;

                /* HiiInvHij <- HiiInv*Hg_ij (=inv(Hg_ii)*Hg_ij)  */
                DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.HiiInvs[this_thread][0], &bsz, &sys.Hg[bsz2*pij],
                      &bsz, &zero, &sys.HiiInvHijs[this_thread][0], &bsz);
                sys.C[pij] = DNRM2(&bsz2, &sys.HiiInvHijs[this_thread][0], &inc);
            }
        }
        
//         #pragma omp for private(i, j, m, ri, rj, nri, nrk, ik, kj, kk, pij, pim, pki, Cik, Ckj, shared)
        for (Int rk=0; rk<nrows; rk++) {
            sys.w[rk] = 0.0;
            kk = rowpts[rk];
            nrk = rowpts[rk+1] - rowpts[rk];
            for (i = 1; i<nrk; i++) {
                pki = kk+i;
                ri = colind[pki];
                nri = rowpts[ri+1] - rowpts[ri];
                for (j=1; j<nri; j++) {
                    pij = rowpts[ri]+j;
                    if (colind[pij] == rk) break;
                }
                ik = pij;
                Cik = sys.C[ik];
                for (j = 1; j<nrk; j++)
                    if (j != i) {
                        kj = kk+j;
                        rj = colind[kj];
                        Ckj = sys.C[kj];

                        shared = 0;
                        for (m=1; m<nri; m++) {     // Determine if j neighbors i
                            pim = rowpts[ri]+m;
                            if (colind[pim] == rj) {
                                shared = 1;
                                break;
                            }
                        }
                        if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                            sys.w[rk] += Cik*Cik*Ckj*Ckj;
                        }
                    }
            }
            sys.w[rk] = sqrt(sys.w[rk]);
        }
//     }
        
        // Check NaN and Inf in sys.w:
        for (i = 0; i < nrows; i++) {
            if (isnan(sys.w[i]) || isinf(sys.w[i]) || (sys.w[i] > bign)) {
                printf("Error in MDF reorder algorithm: sys.w[%d] = %g.\nExecution will be terminated.\n", i, sys.w[i]);
                exit(-1);
            }
        }
        
        for (Int rl=0; rl<nrows; rl++) {
// // //             Int chunkLen, startPt;
// // //             splitArray1D(nrows, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //             pl_local[this_thread] = startPt + indexofSmallestElement(&sys.w[startPt], chunkLen);
            pl = indexofSmallestElement(&sys.w[0], nrows);
// // //             #pragma omp barrier     // TODO: This barrier (or in general the threaded way of computing pl may be too slow)
// // //             #pragma omp single
// // //             {
// // //                 pl = pl_local[0];
// // //                 for (Int i_th = 1; i_th < sys.noThreads; i_th++)
// // //                     if (sys.w[pl_local[i_th]] < sys.w[pl])
// // //                         pl = pl_local[i_th];
                sys.ordered2unordered[rl] = pl;
                sys.w[pl] = bign;
                nrl = rowpts[pl+1] - rowpts[pl];
// // //             }
// // //             #pragma omp for private(i, j, k, m, ri, rj, rk, nri, nrk, pij, pki, plk, pim, ik, kj, Cik, Ckj, shared)
            for (k = 1; k<nrl; k++) {
                plk = rowpts[pl]+k;
                rk = colind[plk];
                if (sys.w[rk] < bign) {
                    sys.w[rk] = 0.0;
                    nrk = rowpts[rk+1] - rowpts[rk];
                    for (i = 1; i<nrk; i++) {
                        pki = rowpts[rk]+i;
                        ri = colind[pki];
                        if (sys.w[ri] < bign) {
                            nri = rowpts[ri+1] - rowpts[ri];
                            for (j = 1; j<nri; j++) {
                                pij = rowpts[ri]+j;
                                if (colind[pij] == rk) break;
                            }
                            ik = pij;
                            Cik = sys.C[ik];
                            for (j = 1; j<nrk; j++) {
                                kj = rowpts[rk]+j;
                                rj = colind[kj];
                                if (j != i && sys.w[rj] < bign) {
                                    Ckj = sys.C[kj];

                                    shared = 0;
                                    for (m = 1; m < nri; m++) {     // Determine if j neighbors i
                                        pim = rowpts[ri] + m;
                                        if (colind[pim] == rj) {
                                            shared = 1;
                                            break;
                                        }
                                    }
                                    if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                                        sys.w[rk] += Cik * Cik * Ckj * Ckj;
                                    }
                                }
                            }
                        }
                    }
                    sys.w[rk] = sqrt(sys.w[rk]);
                }
            }
        }
        
//         #pragma omp parallel for
        for (Int rk=0; rk<nrows; rk++)
            sys.unordered2ordered[sys.ordered2unordered[rk]] = rk;
    
    delete[] pl_local;
}

void computeExactMDFordering(sysstruct &sys, Int preconditioner)
{
    /* TODO: Further optimization: If necessary, this function can be optimized by not recomputing w[i] when needs to be updated (we can just remove the contributions that no longer exist). */

    /* Note:
     * In this function, entity refers to faces in HDG and to face nodes in EDG and IEDG.
     * We assume all blocks are of the same size (i.e., nch (EDG) or nch*npf (UDH) are constant throughout the mesh).
     *      If not, an array indicating the size of each block and an array indicating the starting position each block in the matrix J would be required (easy fix).
     * The complexity of the MDF ordering is O(N^2) due to the indexofSmallestElement operation.
     * */

    char chn = 'N';
    double one = 1.0, zero = 0.0;

    Int nrows, nblks;
    Int *rowpts, *colind;
    if (preconditioner == 0) {
        nrows = sys.numEntities;
        nblks = sys.numBlocks;
        rowpts = &sys.ent2entStart[0];
        colind = &sys.ent2ent[0];
    }
    else if (preconditioner == 1) {
        nrows = sys.BJ_nrows;
        nblks = sys.numBlocks;
        rowpts = &sys.ent2entStart[0];
        colind = &sys.ent2ent[0];
        //nrows = sys.BJ_nrows;
        //nblks = sys.BJ_nblks;
        //rowpts = &sys.BJ_rowpts[0];
        //colind = &sys.BJ_colind[0];
    }
    else
        error("Preconditioner flag not recognized in computeExactMDFordering.\n");
    
    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int inc = 1, i,j,k,m,n, info;
    Int rk, kk, nrk, ri, rj, nri, pij, pki, ik, kj, rl, pl, pim, nrl, plk, shared;
    double bign = 1.0e100, dnr;

    vector<vector<Int> >* ipivs = &sys.ipivs;
    Int lwork = bsz;
    vector<vector<double> >* works = &sys.works;
    Int *pl_local = new Int[sys.noThreads];

//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for private(i, j, m, ri, rj, nri, nrk, ik, kj, kk, pij, pim, pki, dnr, shared, info)
        for (Int rk=0; rk<nrows; rk++) {
            sys.w[rk] = 0.0;
            kk = rowpts[rk];
            nrk = rowpts[rk+1] - rowpts[rk];

            // HkkInv <- inv(Hg_kk)
            DCOPY(&bsz2, &sys.Hg[bsz2*kk], &inc, &sys.HkkInvs[this_thread][0], &inc);
            DGETRF(&bsz, &bsz, &sys.HkkInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &info);
            DGETRI(&bsz, &sys.HkkInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
            for (i = 1; i<nrk; i++) {
                pki = kk+i;
                ri = colind[pki];
                nri = rowpts[ri+1] - rowpts[ri];
                for (j=1; j<nri; j++) {
                    pij = rowpts[ri]+j;
                    if (colind[pij] == rk) break;
                }
                ik = pij;

                /* HikHkkInv <- H_ik*HkkInv (=H_ik*inv(H_kk)) */
                DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.Hg[bsz2*ik], &bsz, &sys.HkkInvs[this_thread][0],
                      &bsz, &zero, &sys.HikHkkInvs[this_thread][0], &bsz);
                for (j = 1; j<nrk; j++)
                    if (j != i) {
                        kj = kk+j;
                        rj = colind[kj];

                        shared = 0;
                        for (m=1; m<nri; m++) {     // Determine if j neighbors i
                            pim = rowpts[ri]+m;
                            if (colind[pim] == rj) {
                                shared = 1;
                                break;
                            }
                        }
                        if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                            /* HikHkkInvHkj <- HikHkkInv*H_kj (=H_ik*inv(H_kk)*H_kj) */
                            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.HikHkkInvs[this_thread][0], &bsz, &sys.Hg[bsz2*kj],
                                  &bsz, &zero, &sys.HikHkkInvHkjs[this_thread][0], &bsz);
                            dnr = DNRM2(&bsz2, &sys.HikHkkInvHkjs[this_thread][0], &inc);
                            sys.w[rk] += dnr*dnr;   // TODO: Is this line correct or are we doing ^4?? (same with algorithm with constraints)
                        }
                    }
            }
            sys.w[rk] = sqrt(sys.w[rk]);
        }
        
        // Check NaN and Inf in sys.w:
        for (i = 0; i < nrows; i++) {
            if (isnan(sys.w[i]) || isinf(sys.w[i]) || (sys.w[i] > bign)) {
                printf("Error in MDF reorder algorithm: sys.w[%d] = %g.\nExecution will be terminated.\n", i, sys.w[i]);
                exit(-1);
            }
        }
        
        for (Int rl=0; rl<nrows; rl++) {
// // //             Int chunkLen, startPt;
// // //             splitArray1D(nrows, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //             pl_local[this_thread] = startPt + indexofSmallestElement(&sys.w[startPt], chunkLen);
            pl = indexofSmallestElement(&sys.w[0], nrows);
// // //             #pragma omp barrier
// // //             #pragma omp single
// // //             {
// // //                 pl = pl_local[0];
// // //                 for (Int i_th = 1; i_th < sys.noThreads; i_th++)
// // //                     if (sys.w[pl_local[i_th]] < sys.w[pl])
// // //                         pl = pl_local[i_th];
                sys.ordered2unordered[rl] = pl;
                sys.w[pl] = bign;
                nrl = rowpts[pl+1] - rowpts[pl];
// // //             }
//             #pragma omp for private(i, j, k, m, kk, ri, rj, rk, nri, nrk, pij, pki, plk, pim, ik, kj, dnr, shared, info)
            for (k = 1; k<nrl; k++) {
                plk = rowpts[pl]+k;
                rk = colind[plk];
                kk = rowpts[rk];
                if (sys.w[rk] < bign) {
                    sys.w[rk] = 0.0;
                    nrk = rowpts[rk+1] - rowpts[rk];

                    // HkkInv <- inv(Hg_kk)
                    DCOPY(&bsz2, &sys.Hg[bsz2*kk], &inc, &sys.HkkInvs[this_thread][0], &inc);
                    DGETRF(&bsz, &bsz, &sys.HkkInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &info);
                    DGETRI(&bsz, &sys.HkkInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
                    for (i = 1; i<nrk; i++) {
                        pki = rowpts[rk]+i;
                        ri = colind[pki];
                        if (sys.w[ri] < bign) {
                            nri = rowpts[ri+1] - rowpts[ri];
                            for (j = 1; j< nri; j++) {
                                pij = rowpts[ri]+j;
                                if (colind[pij] == rk) break;
                            }
                            ik = pij;
                            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.Hg[bsz2*ik], &bsz, &sys.HkkInvs[this_thread][0],
                                  &bsz, &zero, &sys.HikHkkInvs[this_thread][0], &bsz);
                            for (j = 1; j<nrk; j++) {
                                kj = rowpts[rk]+j;
                                rj = colind[kj];
                                if (j != i && sys.w[rj] < bign) {
                                    shared = 0;
                                    for (m = 1; m < nri; m++) {     // Determine if j neighbors i
                                        pim = rowpts[ri] + m;
                                        if (colind[pim] == rj) {
                                            shared = 1;
                                            break;
                                        }
                                    }
                                    if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                                        DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.HikHkkInvs[this_thread][0], &bsz, &sys.Hg[bsz2*kj],
                                              &bsz, &zero, &sys.HikHkkInvHkjs[this_thread][0], &bsz);
                                        dnr = DNRM2(&bsz2, &sys.HikHkkInvHkjs[this_thread][0], &inc);
                                        sys.w[rk] += dnr*dnr;
                                    }
                                }
                            }
                        }
                    }
                    sys.w[rk] = sqrt(sys.w[rk]);
                }
            }
        }

//         #pragma omp for
        for (Int rk=0; rk<nrows; rk++)
            sys.unordered2ordered[sys.ordered2unordered[rk]] = rk;
//     }
    
    delete[] pl_local;
}

void computeApproximateMDForderingWithConstraints(sysstruct &sys, Int preconditioner)
{
    /* TODO: Further optimization: If necessary, this function can be optimized by not recomputing w[i] when needs to be updated (we can just remove the contributions that no longer exist). */

    // Constraints: Neighboring entities are ordered last

    /* Note:
     * In this function, entity refers to faces in HDG and to face nodes in EDG and IEDG.
     * We assume all blocks are of the same size (i.e., nch (EDG) or nch*npf (UDH) are constant throughout the mesh).
     *      If not, an array indicating the size of each block and an array indicating the starting position each block in the matrix J would be required (easy fix).
     * The complexity of the MDF ordering is O(N^2) due to the indexofSmallestElement operation.
     * */

    if (sys.nproc == 1) {
        error("\nApproximate MDF ordering with constraints only valid for parallel RAS preconditioner.\n");
    }
    else if (preconditioner != 0) {
        error("\nApproximate MDF ordering with constraints only valid for parallel RAS preconditioner.\n");
    }

    char chn = 'N';
    double one = 1.0, zero = 0.0;

    Int nrows, nrowsBJ, nblks;
    Int *rowpts, *colind;
    nrows = sys.numEntities;
    nblks = sys.numBlocks;
    rowpts = &sys.ent2entStart[0];
    colind = &sys.ent2ent[0];
    nrowsBJ = sys.BJ_nrows;

    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int inc = 1, i,j,k,m,n, info;
    Int rk, kk, nrk, ri, rj, nri, pij, pki, ii, ik, kj, rl, pl, pim, nrl, plk, shared;
    double Cik, Dik, Ckj, bign = 1.0e100;

    vector<vector<Int> >* ipivs = &sys.ipivs;
    Int lwork = bsz;
    vector<vector<double> >* works = &sys.works;
    Int *pl_local = new Int[sys.noThreads];

//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for private(info, j, ii, nri, pij)
        for (Int ri=0; ri<nrows; ri++) {
            ii = rowpts[ri];
            nri = rowpts[ri+1] - rowpts[ri];

            //  HiiInv <- inv(Hg_ii)
            DCOPY(&bsz2, &sys.Hg[bsz2*ii], &inc, &sys.HiiInvs[this_thread][0], &inc);
            DGETRF(&bsz, &bsz, &sys.HiiInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &info);
            DGETRI(&bsz, &sys.HiiInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
            for (j=0; j<nri; j++) {
                pij = ii+j;

                // HiiInvHij <- HiiInv*Hg_ij (=inv(Hg_ii)*Hg_ij)
                DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.HiiInvs[this_thread][0], &bsz, &sys.Hg[bsz2*pij],
                      &bsz, &zero, &sys.HiiInvHijs[this_thread][0], &bsz);
                sys.C[pij] = DNRM2(&bsz2, &sys.HiiInvHijs[this_thread][0], &inc);
            }
        }

//         #pragma omp for private(i, j, m, ri, rj, nri, nrk, ik, kj, kk, pij, pim, pki, Cik, Ckj, shared)
        for (Int rk=0; rk<nrows; rk++) {
            sys.w[rk] = 0.0;
            kk = rowpts[rk];
            nrk = rowpts[rk+1] - rowpts[rk];
            for (i = 1; i<nrk; i++) {
                pki = kk+i;
                ri = colind[pki];
                nri = rowpts[ri+1] - rowpts[ri];
                for (j=1; j<nri; j++) {
                    pij = rowpts[ri]+j;
                    if (colind[pij] == rk) break;
                }
                ik = pij;
                Cik = sys.C[ik];
                for (j = 1; j<nrk; j++)
                    if (j != i) {
                        kj = kk+j;
                        rj = colind[kj];
                        Ckj = sys.C[kj];

                        shared = 0;
                        for (m=1; m<nri; m++) {     // Determine if j neighbors i
                            pim = rowpts[ri]+m;
                            if (colind[pim] == rj) {
                                shared = 1;
                                break;
                            }
                        }
                        if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                            sys.w[rk] += Cik*Cik*Ckj*Ckj;
                        }
                    }
            }
            sys.w[rk] = sqrt(sys.w[rk]);
        }
//     }
        
        // Check NaN and Inf in sys.w:
        for (i = 0; i < nrows; i++) {
            if (isnan(sys.w[i]) || isinf(sys.w[i]) || (sys.w[i] > bign)) {
                printf("Error in MDF reorder algorithm: sys.w[%d] = %g.\nExecution will be terminated.\n", i, sys.w[i]);
                exit(-1);
            }
        }
        
        for (Int rl=0; rl<nrows; rl++) {
// // //             Int chunkLen, startPt;
// // //             if (rl < nrowsBJ) {
// // //                 splitArray1D(nrowsBJ, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //                 pl_local[this_thread] = startPt + indexofSmallestElement(&sys.w[startPt], chunkLen);
// // //             }
// // //             else {
// // //                 splitArray1D(nrows-nrowsBJ, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //                 pl_local[this_thread] = nrowsBJ+startPt + indexofSmallestElement(&sys.w[nrowsBJ+startPt], chunkLen);
// // //             }
            if (rl < nrowsBJ)
                pl = indexofSmallestElement(&sys.w[0], nrowsBJ);
            else
                pl = nrowsBJ + indexofSmallestElement(&sys.w[nrowsBJ], nrows-nrowsBJ);
// // //             #pragma omp barrier
// // //             #pragma omp single
// // //             {
// // //                 pl = pl_local[0];
// // //                 for (Int i_th = 1; i_th < sys.noThreads; i_th++)
// // //                     if (sys.w[pl_local[i_th]] < sys.w[pl])
// // //                         pl = pl_local[i_th];
                sys.ordered2unordered[rl] = pl;
                sys.w[pl] = bign;
                nrl = rowpts[pl+1] - rowpts[pl];
// // //             }
// // //             #pragma omp for private(i, j, k, m, ri, rj, rk, nri, nrk, pij, pki, plk, pim, ik, kj, Cik, Ckj, shared)
            for (k = 1; k<nrl; k++) {
                plk = rowpts[pl]+k;
                rk = colind[plk];
                if (sys.w[rk] < bign) {
                    sys.w[rk] = 0.0;
                    nrk = rowpts[rk+1] - rowpts[rk];
                    for (i = 1; i<nrk; i++) {
                        pki = rowpts[rk]+i;
                        ri = colind[pki];
                        if (sys.w[ri] < bign) {
                            nri = rowpts[ri+1] - rowpts[ri];
                            for (j = 1; j<nri; j++) {
                                pij = rowpts[ri]+j;
                                if (colind[pij] == rk) break;
                            }
                            ik = pij;
                            Cik = sys.C[ik];
                            for (j = 1; j<nrk; j++) {
                                kj = rowpts[rk]+j;
                                rj = colind[kj];
                                if (j != i && sys.w[rj] < bign) {
                                    Ckj = sys.C[kj];

                                    shared = 0;
                                    for (m = 1; m < nri; m++) {     // Determine if j neighbors i
                                        pim = rowpts[ri] + m;
                                        if (colind[pim] == rj) {
                                            shared = 1;
                                            break;
                                        }
                                    }
                                    if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                                        sys.w[rk] += Cik * Cik * Ckj * Ckj;
                                    }
                                }
                            }
                        }
                    }
                    sys.w[rk] = sqrt(sys.w[rk]);
                }
            }
        }

//         #pragma omp for
        for (Int rk=0; rk<nrows; rk++)
            sys.unordered2ordered[sys.ordered2unordered[rk]] = rk;
    
    delete[] pl_local;
}

void computeExactMDForderingWithConstraints(sysstruct &sys, Int preconditioner)
{
    /* TODO: Further optimization: If necessary, this function can be optimized by not recomputing w[i] when needs to be updated (we can just remove the contributions that no longer exist). */

    // Constraints: Neighboring entities are ordered last

     /* Note:
     * In this function, entity refers to faces in HDG and to face nodes in EDG and IEDG.
     * We assume all blocks are of the same size (i.e., nch (EDG) or nch*npf (UDH) are constant throughout the mesh).
     *      If not, an array indicating the size of each block and an array indicating the starting position each block in the matrix J would be required (easy fix).
     * The complexity of the MDF ordering is O(N^2) due to the indexofSmallestElement operation.
     * */
    
    if (sys.nproc == 1) {
        error("\nExact MDF ordering with constraints only valid for parallel RAS preconditioner.\n");
    }
    else if (preconditioner != 0) {
        error("\nExact MDF ordering with constraints only valid for parallel RAS preconditioner.\n");
    }

    char chn = 'N';
    double one = 1.0, zero = 0.0;

    Int nrows, nrowsBJ, nblks;
    Int *rowpts, *colind;
    nrows = sys.numEntities;
    nblks = sys.numBlocks;
    rowpts = &sys.ent2entStart[0];
    colind = &sys.ent2ent[0];
    nrowsBJ = sys.BJ_nrows;

    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int inc = 1, i,j,k,m,n, info;
    Int rk, kk, nrk, ri, rj, nri, pij, pki, ik, kj, rl, pl, pim, nrl, plk, shared;
    double bign = 1.0e100, dnr;

    Int *ipiv = &sys.ipiv[0];
    vector<vector<Int> >* ipivs = &sys.ipivs;
    Int lwork = bsz;
    double *work = &sys.work[0];
    vector<vector<double> >* works = &sys.works;
    Int *pl_local = new Int[sys.noThreads];

//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for private(i, j, m, ri, rj, nri, nrk, ik, kj, kk, pij, pim, pki, dnr, shared, info)
        for (Int rk=0; rk<nrows; rk++) {
            sys.w[rk] = 0.0;
            kk = rowpts[rk];
            nrk = rowpts[rk+1] - rowpts[rk];

            // HkkInv <- inv(Hg_kk)
            DCOPY(&bsz2, &sys.Hg[bsz2*kk], &inc, &sys.HkkInvs[this_thread][0], &inc);
            DGETRF(&bsz, &bsz, &sys.HkkInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &info);
            DGETRI(&bsz, &sys.HkkInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
            for (i = 1; i<nrk; i++) {
                pki = kk+i;
                ri = colind[pki];
                nri = rowpts[ri+1] - rowpts[ri];
                for (j=1; j<nri; j++) {
                    pij = rowpts[ri]+j;
                    if (colind[pij] == rk) break;
                }
                ik = pij;

                // HikHkkInv <- H_ik*HkkInv (=H_ik*inv(H_kk))
                DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.Hg[bsz2*ik], &bsz, &sys.HkkInvs[this_thread][0],
                      &bsz, &zero, &sys.HikHkkInvs[this_thread][0], &bsz);
                for (j = 1; j<nrk; j++)
                    if (j != i) {
                        kj = kk+j;
                        rj = colind[kj];

                        shared = 0;
                        for (m=1; m<nri; m++) {     // Determine if j neighbors i
                            pim = rowpts[ri]+m;
                            if (colind[pim] == rj) {
                                shared = 1;
                                break;
                            }
                        }
                        if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                            /* HikHkkInvHkj <- HikHkkInv*H_kj (=H_ik*inv(H_kk)*H_kj) */
                            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.HikHkkInvs[this_thread][0], &bsz, &sys.Hg[bsz2*kj],
                                  &bsz, &zero, &sys.HikHkkInvHkjs[this_thread][0], &bsz);
                            dnr = DNRM2(&bsz2, &sys.HikHkkInvHkjs[this_thread][0], &inc);
                            sys.w[rk] += dnr*dnr;
                        }
                    }
            }
            sys.w[rk] = sqrt(sys.w[rk]);
        }

        // Check NaN and Inf in sys.w:
        for (i = 0; i < nrows; i++) {
            if (isnan(sys.w[i]) || isinf(sys.w[i]) || (sys.w[i] > bign)) {
                printf("Error in MDF reorder algorithm: sys.w[%d] = %g.\nExecution will be terminated.\n", i, sys.w[i]);
                exit(-1);
            }
        }
        
        for (Int rl=0; rl<nrows; rl++) {
// // //             Int chunkLen, startPt;
// // //             if (rl < nrowsBJ) {
// // //                 splitArray1D(nrowsBJ, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //                 pl_local[this_thread] = startPt + indexofSmallestElement(&sys.w[startPt], chunkLen);
// // //             }
// // //             else {
// // //                 splitArray1D(nrows-nrowsBJ, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //                 pl_local[this_thread] = nrowsBJ+startPt + indexofSmallestElement(&sys.w[nrowsBJ+startPt], chunkLen);
// // //             }
            if (rl < nrowsBJ)
                pl = indexofSmallestElement(&sys.w[0], nrowsBJ);
            else
                pl = nrowsBJ + indexofSmallestElement(&sys.w[nrowsBJ], nrows-nrowsBJ);
// // //             #pragma omp barrier
// // //             #pragma omp single
// // //             {
// // //                 pl = pl_local[0];
// // //                 for (Int i_th = 1; i_th < sys.noThreads; i_th++)
// // //                     if (sys.w[pl_local[i_th]] < sys.w[pl])
// // //                         pl = pl_local[i_th];
                sys.ordered2unordered[rl] = pl;
                sys.w[pl] = bign;
                nrl = rowpts[pl+1] - rowpts[pl];
// // //             }
//             #pragma omp for private(i, j, k, m, kk, ri, rj, rk, nri, nrk, pij, pki, plk, pim, ik, kj, dnr, shared, info)
            for (k = 1; k<nrl; k++) {
                plk = rowpts[pl]+k;
                rk = colind[plk];
                kk = rowpts[rk];
                if (sys.w[rk] < bign) {
                    sys.w[rk] = 0.0;
                    nrk = rowpts[rk+1] - rowpts[rk];

                    // HkkInv <- inv(Hg_kk)
                    DCOPY(&bsz2, &sys.Hg[bsz2*kk], &inc, &sys.HkkInvs[this_thread][0], &inc);
                    DGETRF(&bsz, &bsz, &sys.HkkInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &info);
                    DGETRI(&bsz, &sys.HkkInvs[this_thread][0], &bsz, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
                    for (i = 1; i<nrk; i++) {
                        pki = rowpts[rk]+i;
                        ri = colind[pki];
                        if (sys.w[ri] < bign) {
                            nri = rowpts[ri+1] - rowpts[ri];
                            for (j = 1; j< nri; j++) {
                                pij = rowpts[ri]+j;
                                if (colind[pij] == rk) break;
                            }
                            ik = pij;
                            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.Hg[bsz2*ik], &bsz, &sys.HkkInvs[this_thread][0],
                                  &bsz, &zero, &sys.HikHkkInvs[this_thread][0], &bsz);
                            for (j = 1; j<nrk; j++) {
                                kj = rowpts[rk]+j;
                                rj = colind[kj];
                                if (j != i && sys.w[rj] < bign) {
                                    shared = 0;
                                    for (m = 1; m < nri; m++) {     // Determine if j neighbors i
                                        pim = rowpts[ri] + m;
                                        if (colind[pim] == rj) {
                                            shared = 1;
                                            break;
                                        }
                                    }
                                    if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                                        DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.HikHkkInvs[this_thread][0], &bsz, &sys.Hg[bsz2*kj],
                                              &bsz, &zero, &sys.HikHkkInvHkjs[this_thread][0], &bsz);
                                        dnr = DNRM2(&bsz2, &sys.HikHkkInvHkjs[this_thread][0], &inc);
                                        sys.w[rk] += dnr*dnr;
                                    }
                                }
                            }
                        }
                    }
                    sys.w[rk] = sqrt(sys.w[rk]);
                }
            }
        }

//         #pragma omp for
        for (Int rk=0; rk<nrows; rk++)
            sys.unordered2ordered[sys.ordered2unordered[rk]] = rk;
//     }
    
    delete[] pl_local;
}

void postprocessOrderingForBILU0(sysstruct &sys, Int preconditioner)
{
    // TODO: This function is hard to thread and the benefit may be small.
    
    Int nrows, nblks;
    Int *rowpts, *colind;
    if (preconditioner == 1) {
        nrows = sys.BJ_nrows;
        nblks = sys.numBlocks;
        rowpts = &sys.ent2entStart[0];
        colind = &sys.ent2ent[0];
        //nrows = sys.BJ_nrows;
        //nblks = sys.BJ_nblks;
        //rowpts = &sys.BJ_rowpts[0];
        //colind = &sys.BJ_colind[0];
    }
    else
        error("Preconditioner flag not recognized in postprocessOrderingForBILU0.\n");
    
    Int i, j, i_l, i_u;
    Int rj, rjj, ri, rii, pj, pjj, pj_l, pj_u, pjji, nrj;
    
    // Compute sys.oent2entStart and sys.oent2ent
    sys.oent2entStart[0] = 0;
    for (rj=0; rj<nrows; rj++) {
        rjj = sys.ordered2unordered[rj];
        pj = sys.oent2entStart[rj];
        pjj = rowpts[rjj];
        nrj = rowpts[rjj+1]-rowpts[rjj];

        sys.oent2entStart[rj+1] = sys.oent2entStart[rj] + nrj;
        for (i=0; i<nrj; i++)
            sys.oent2ent[pj+i] = colind[pjj+i];
    }

    // Compute sys.LUoent2ent (part associated to L) and sys.Loent2entStart
    sys.Loent2entStart[0] = 0;
    for (rj=0; rj<nrows; rj++) {
        rjj = sys.ordered2unordered[rj];

        pj_l = sys.Loent2entStart[rj];
        pjj = rowpts[rjj];

        nrj = rowpts[rjj+1]-rowpts[rjj];

        i_l = 0;
        for (i=1; i<nrj; i++) {
            pjji = pjj+i;
            rii = colind[pjji];
            ri = sys.unordered2ordered[rii];
            if (ri < rj) {
                sys.LUoent2ent[pj_l+i_l] = rii;
                i_l ++;
            }
        }
        sys.Loent2entStart[rj+1] = sys.Loent2entStart[rj] + i_l;
    }

    // Compute sys.Uoent2entStart
    sys.Uoent2entStart[nrows] = sys.Loent2entStart[nrows];
    for (j=0; j<nrows; j++) {
        rj = nrows-j-1;
        rjj = sys.ordered2unordered[rj];

        pjj = rowpts[rjj];

        nrj = rowpts[rjj+1]-rowpts[rjj];

        if (nrj > 0) {
            i_u = 1;        // Start at 1 to consider the diagonal block
        }
        else {              // This is a ghost (missing) entity in the mesh
            printf("\nWARNING: Ghost entities have been detected in the mesh.\n\n");
            i_u = 0;
        }

        for (i=1; i<nrj; i++) {
            pjji = pjj+i;
            rii = colind[pjji];
            ri = sys.unordered2ordered[rii];
            if (ri > rj) {
                i_u ++;
            }
        }
        sys.Uoent2entStart[rj] = sys.Uoent2entStart[rj+1] + i_u;
    }

    // Compute sys.LUoent2ent (part associated to U)
    for (j=0; j<nrows; j++) {
        rj = nrows-j-1;
        rjj = sys.ordered2unordered[rj];

        pj_u = sys.Uoent2entStart[rj+1];
        pjj = rowpts[rjj];

        nrj = rowpts[rjj+1] - rowpts[rjj];

        i_u = 0;
        for (i=1; i<nrj; i++) {
            pjji = pjj+i;
            rii = colind[pjji];
            ri = sys.unordered2ordered[rii];
            if (ri > rj) {
                sys.LUoent2ent[pj_u+i_u] = rii;
                i_u ++;
            }
        }
        sys.LUoent2ent[pj_u+i_u] = rjj;
    }
}

void postprocessOrderingForBILUk(sysstruct &sys)
{
    // TODO: This function is hard to thread and the benefit may be small.
    
    Int nrowsP = sys.nrowsP;
    Int nblksP = sys.nblksP;
    Int *rowptsP = &sys.rowptsP[0];
    Int *colindP = &sys.colindP[0];
    
    Int i, j, i_l, i_u;
    Int rj, rjj, ri, rii, pj, pjj, pj_l, pj_u, pjji, nrj;
    
    // Compute sys.oent2entStart and sys.oent2ent
    sys.oent2entStart[0] = 0;
    for (rj=0; rj<nrowsP; rj++) {
        rjj = sys.ordered2unordered[rj];
        pj = sys.oent2entStart[rj];
        pjj = rowptsP[rjj];
        nrj = rowptsP[rjj+1]-rowptsP[rjj];

        sys.oent2entStart[rj+1] = sys.oent2entStart[rj] + nrj;
        for (i=0; i<nrj; i++)
            sys.oent2ent[pj+i] = colindP[pjj+i];
    }

    // Compute sys.LUoent2ent (part associated to L) and sys.Loent2entStart
    sys.Loent2entStart[0] = 0;
    for (rj=0; rj<nrowsP; rj++) {
        rjj = sys.ordered2unordered[rj];
        
        pj_l = sys.Loent2entStart[rj];
        pjj = rowptsP[rjj];
        
        nrj = rowptsP[rjj+1]-rowptsP[rjj];
        
        i_l = 0;
        for (i=1; i<nrj; i++) {
            pjji = pjj+i;
            rii = colindP[pjji];
            ri = sys.unordered2ordered[rii];
            if (ri < rj) {
                sys.LUoent2ent[pj_l+i_l] = rii;
                i_l ++;
            }
        }
        sys.Loent2entStart[rj+1] = sys.Loent2entStart[rj] + i_l;
    }

    // Compute sys.Uoent2entStart
    sys.Uoent2entStart[nrowsP] = sys.Loent2entStart[nrowsP];
    for (j=0; j<nrowsP; j++) {
        rj = nrowsP-j-1;
        rjj = sys.ordered2unordered[rj];

        pjj = rowptsP[rjj];

        nrj = rowptsP[rjj+1]-rowptsP[rjj];

        if (nrj > 0) {
            i_u = 1;        // Start at 1 to consider the diagonal block
        }
        else {              // This is a ghost (missing) entity in the mesh
            printf("\nWARNING: Ghost entities have been detected in the mesh.\n\n");
            i_u = 0;
        }

        for (i=1; i<nrj; i++) {
            pjji = pjj+i;
            rii = colindP[pjji];
            ri = sys.unordered2ordered[rii];
            if (ri > rj) {
                i_u ++;
            }
        }
        sys.Uoent2entStart[rj] = sys.Uoent2entStart[rj+1] + i_u;
    }
    
    // Compute sys.LUoent2ent (part associated to U)
    for (j=0; j<nrowsP; j++) {
        rj = nrowsP-j-1;
        rjj = sys.ordered2unordered[rj];
        
        pj_u = sys.Uoent2entStart[rj+1];
        pjj = rowptsP[rjj];
        
        nrj = rowptsP[rjj+1] - rowptsP[rjj];

        i_u = 0;
        for (i=1; i<nrj; i++) {
            pjji = pjj+i;
            rii = colindP[pjji];
            ri = sys.unordered2ordered[rii];
            if (ri > rj) {
                sys.LUoent2ent[pj_u+i_u] = rii;
                i_u ++;
            }
        }
        sys.LUoent2ent[pj_u+i_u] = rjj;
    }
}

void copyHg2MgWithMgFormat(sysstruct &sys)
{
    Int inc = 1;
    Int nrowsP = sys.nrowsP;
    
    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    
    Int *rowptsJ = &sys.ent2entStart[0];
    Int *rowptsP = &sys.rowptsP[0];
    Int *colindP = &sys.colindP[0];
    
    for (Int rjj=0; rjj<nrowsP; rjj++) {
        Int rj = sys.unordered2ordered[rjj];
        
        Int pjP   = sys.Uoent2entStart[rj]-1;
        Int pjP_l = sys.Loent2entStart[rj];
        Int pjP_u = sys.Uoent2entStart[rj+1];
        
        Int pjjJ  = rowptsJ[rjj];
        Int pjjP  = rowptsP[rjj];
        
        Int nrjJ  = rowptsJ[rjj+1]-rowptsJ[rjj];
        Int nrjP  = rowptsP[rjj+1]-rowptsP[rjj];
        Int nrj_l = sys.Loent2entStart[rj+1]-sys.Loent2entStart[rj];        // Neighbors in L
        Int nrj_u = sys.Uoent2entStart[rj]-sys.Uoent2entStart[rj+1];        // Neighbors in U (includes the diagonal block)
        
        // Copy diagonal block:
        DCOPY(&bsz2, &sys.Hg[bsz2*pjjJ], &inc, &sys.Mg[bsz2*pjP], &inc);
        
        // Copy off-diagonal blocks:
        Int i_u = 0, i_l = 0;
        Int pjiP_l = pjP_l;
        Int pjiP_u = pjP_u;
        Int pjjiJ = pjjJ + 1;
        for (Int i=1; i<nrjP; i++) {
            Int rii = colindP[pjjP+i];
            Int ri = sys.unordered2ordered[rii];
            
            if (ri > rj) {
                // U block
                i_u ++;
                if (! sys.MgMinusHg[rowptsP[rjj]+i])
                    std::fill(&sys.Mg[bsz2*pjiP_u], &sys.Mg[bsz2*(pjiP_u+1)], 0.0);
                else {
                    DCOPY(&bsz2, &sys.Hg[bsz2*pjjiJ], &inc, &sys.Mg[bsz2*pjiP_u], &inc);
                    pjjiJ += 1;
                }
                pjiP_u += 1;
            }
            else if (ri < rj) {
                // L block
                i_l ++;
                if (! sys.MgMinusHg[rowptsP[rjj]+i])
                    std::fill(&sys.Mg[bsz2*pjiP_l], &sys.Mg[bsz2*(pjiP_l+1)], 0.0);
                else {
                    DCOPY(&bsz2, &sys.Hg[bsz2*pjjiJ], &inc, &sys.Mg[bsz2*pjiP_l], &inc);
                    pjjiJ += 1;
                }
                pjiP_l += 1;
            }
            else {
                printf("%d,%d,%d,%d,%d\n",ri,rii,rj,rjj,i);
                error("Error No. 1 in copyHg2MgWithMgFormat.\n");
            }
        }
        
        if ((i_l != nrj_l || i_u != (nrj_u-1)) && nrjP > 0) {
            printf("Error No. 2 in copyHg2MgWithMgFormat. rjj= %d, rj = %d, i_l = %d, nrj_l = %d, i_u = %d, nrj_u = %d, my_rank = %d\n.", rjj, rj, i_l, nrj_l, i_u, nrj_u, sys.my_rank);
            error("\n");
        }
    }
}

void computeBILU0(sysstruct &sys, Int preconditioner)
{
    // Mg is stored in memory in the same order as is required for the BILU0 solve (L first, then U)

    Int nrows, nblks;
    Int *rowpts, *colind;
    if (preconditioner == 1) {
        nrows = sys.BJ_nrows;
        nblks = sys.numBlocks;
        rowpts = &sys.ent2entStart[0];
        colind = &sys.ent2ent[0];
        //nrows = sys.BJ_nrows;
        //nblks = sys.BJ_nblks;
        //rowpts = &sys.BJ_rowpts[0];
        //colind = &sys.BJ_colind[0];
    }
    else
        error("Preconditioner flag not recognized in computeBILU0.\n");
    
    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int inc = 1, i,j,k,l,m,n, info, lenCopy;

    Int rj, rjj, ri, rii, rl, rll;
    Int pj, pj_l, pj_u, pjj, pi, pi_l, pi_u, pii, pjji, pji, pji_l, pji_u, pij, pij_l, pij_u, pil, pil_l, pil_u, piil, pjjm, pjm, pjm_u;
    Int nrj, nrj_l, nrj_u, nri, nri_l, nri_u;
    Int i_l, i_u;

    char chn = 'N';
    double one = 1.0, zero = 0.0, minusone = -1.0;

    Int *ipiv = &sys.ipiv[0];
    Int lwork = bsz;
    double *work = &sys.work[0];
    double *Lij = &sys.HiiInv[0];

    // Copy Hg into Mg with the desired storage format
    for (rjj=0; rjj<nrows; rjj++) {
        rj = sys.unordered2ordered[rjj];

        pj = sys.Uoent2entStart[rj]-1;
        pj_l = sys.Loent2entStart[rj];
        pj_u = sys.Uoent2entStart[rj+1];
        pjj = rowpts[rjj];

        nrj = rowpts[rjj+1]-rowpts[rjj];
        nrj_l = sys.Loent2entStart[rj+1]-sys.Loent2entStart[rj];        // Neighbors in L
        nrj_u = sys.Uoent2entStart[rj]-sys.Uoent2entStart[rj+1];        // Neighbors in U (includes the diagonal block)

        // Copy diagonal block
        DCOPY(&bsz2, &sys.Hg[bsz2*pjj], &inc, &sys.Mg[bsz2*pj], &inc);

        // Copy non-diagonal blocks
        i_u = 0;
        i_l = 0;
        for (i=1; i<nrj; i++) {
            pjji = pjj+i;
            rii = colind[pjji];
            ri = sys.unordered2ordered[rii];

            if (ri > rj) {
                // U block
                pji_u = pj_u + i_u;
                i_u ++;
                DCOPY(&bsz2, &sys.Hg[bsz2*pjji], &inc, &sys.Mg[bsz2*pji_u], &inc);
            }
            else if (ri < rj) {
                // L block
                pji_l = pj_l + i_l;
                i_l ++;
                DCOPY(&bsz2, &sys.Hg[bsz2*pjji], &inc, &sys.Mg[bsz2*pji_l], &inc);
            }
            else {
                error("Error No. 1 in computeBILU0\n");
            }
        }

        if ((i_l != nrj_l || i_u != (nrj_u-1)) && nrj > 0) {
            printf("Error No. 2 in computeBILU0. rjj= %d, rj = %d, i_l = %d, nrj_l = %d, i_u = %d, nrj_u = %d, my_rank = %d\n.", rjj, rj, i_l, nrj_l, i_u, nrj_u, sys.my_rank);
            error("\n");
        }
    }

    // Compute BILU0 preconditioner
    for (rj=0; rj<nrows; rj++) {
        rjj = sys.ordered2unordered[rj];
        pj = sys.Uoent2entStart[rj]-1;
        pj_l = sys.Loent2entStart[rj];
        pj_u = sys.Uoent2entStart[rj+1];
        pjj = rowpts[rjj];

        // Compute inverse of block diagonal of Hg using LU decomposition
        DGETRF(&bsz,&bsz,&sys.Mg[bsz2*pj],&bsz,ipiv,&info);
        DGETRI(&bsz,&sys.Mg[bsz2*pj],&bsz,ipiv,work,&lwork,&info);

        nrj = rowpts[rjj+1]-rowpts[rjj];
        nrj_l = sys.Loent2entStart[rj+1]-sys.Loent2entStart[rj];        // Neighbors in L
        nrj_u = sys.Uoent2entStart[rj]-sys.Uoent2entStart[rj+1];        // Neighbors in U (including the diagonal block)
        for (i=0; i<(nrj_u-1); i++) {
            pji_u = pj_u+i;
            rii = sys.LUoent2ent[pji_u];
            ri = sys.unordered2ordered[rii];

            pi = sys.Uoent2entStart[ri]-1;
            pi_l = sys.Loent2entStart[ri];
            pi_u = sys.Uoent2entStart[ri+1];

            nri = rowpts[rii+1]-rowpts[rii];
            nri_l = sys.Loent2entStart[ri+1]-sys.Loent2entStart[ri];
            nri_u = sys.Uoent2entStart[ri]-sys.Uoent2entStart[ri+1];
            for (l=0; l<nri_l; l++) {
                if (sys.LUoent2ent[pi_l+l] == rjj) break;
            }
            if (l==nri_l && l!=0) {
                error("Error No. 3 in computeBILU0\n");
            }
            pij_l = pi_l + l;

            // Lij <- Uij*inv(Ujj)
            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.Mg[bsz2*pij_l], &bsz, &sys.Mg[bsz2*pj],
                  &bsz, &zero, Lij, &bsz);
            DCOPY(&bsz2, Lij, &inc, &sys.Mg[bsz2*pij_l], &inc);

            // Uii <- Uii - Lij*Uji
            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pij_l], &bsz, &sys.Mg[bsz2*pji_u],
                  &bsz, &one, &sys.Mg[bsz2*pi], &bsz);

            // L_il <- L_il - L_ij*U_jm (m "=" l)
            for (l=0; l<nri_l; l++) {
                pil_l = pi_l+l;
                rll  = sys.LUoent2ent[pil_l];
                rl = sys.unordered2ordered[rll];

                if (rl > rj) {
                    for (m=0; m<(nrj_u-1); m++) {
                        if (rll == sys.LUoent2ent[pj_u+m]) {
                            pjm_u = pj_u+m;
                            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pij_l], &bsz, &sys.Mg[bsz2*pjm_u],
                                  &bsz, &one, &sys.Mg[bsz2*pil_l], &bsz);
                        }
                    }
                }
            }

            // U_il <- U_il - L_ij*U_jm (m "=" l)
            for (l=0; l<(nri_u-1); l++) {
                pil_u = pi_u+l;
                rll  = sys.LUoent2ent[pil_u];
                rl = sys.unordered2ordered[rll];
                
                if (rl > rj) {
                    for (m=0; m<(nrj_u-1); m++) {
                        if (rll == sys.LUoent2ent[pj_u+m]) {
                            pjm_u = pj_u+m;
                            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pij_l], &bsz, &sys.Mg[bsz2*pjm_u],
                                  &bsz, &one, &sys.Mg[bsz2*pil_u], &bsz);
                        }
                    }
                }
            }
        }
    }
}

void computeBILUk(sysstruct &sys)
{
    Int inc = 1, info;
    double one = 1.0, zero = 0.0, minusone = -1.0;
    char chn = 'N';
    
    Int nrowsP = sys.nrowsP;
    Int nblksP = sys.nblksP;
    Int *rowptsP = &sys.rowptsP[0];
    Int *colindP = &sys.colindP[0];
    
    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    
    Int *ipiv = &sys.ipiv[0];
    double *work = &sys.work[0];
    Int lwork = bsz;
    double *Lij = &sys.HiiInv[0];
    
    // Copy Hg into Mg with the desired storage format
    copyHg2MgWithMgFormat(sys);
    
    // Improve diagonal dominance if necessary:
    double epsDiagDom = 0.0;
    double minInftyCondNumber = 1.0e30;
    if (epsDiagDom > 0 || minInftyCondNumber < 1.0e9) {
        for (Int rj=0; rj<nrowsP; rj++) {
            int this_thread = 0;
            
            Int pj = sys.Uoent2entStart[rj]-1;
            Int pj_l = sys.Loent2entStart[rj];
            Int pj_u = sys.Uoent2entStart[rj+1];
            
            Int nrj_l = sys.Loent2entStart[rj+1]-sys.Loent2entStart[rj];        // Neighbors in L
            Int nrj_u = sys.Uoent2entStart[rj]-sys.Uoent2entStart[rj+1];    
            
            // For block algorithm, we use ||A_{ii}^{-1}||^{-1} instead of ||A_{ii}||
            double magDiagBlk = 0.0;
            DCOPY(&bsz2, &sys.Mg[bsz2*pj], &inc, &sys.HiiInvs[this_thread][0], &inc);
            DGETRF(&bsz,&bsz,&sys.HiiInvs[this_thread][0],&bsz,ipiv,&info);
            DGETRI(&bsz,&sys.HiiInvs[this_thread][0],&bsz,ipiv,work,&lwork,&info);
            for (Int j=0; j<bsz2; j++)
                magDiagBlk += abs(sys.HiiInvs[this_thread][j]);
            magDiagBlk = 1.0 / magDiagBlk;
            
            double magPureDiag = 0.0;
            for (Int j=0; j<bsz; j++)
                magPureDiag += abs(sys.Mg[bsz2*pj+j*(bsz+1)]);
            
            double magOffDiagBlk = 0.0;
            for (Int j=0; j<(nrj_u-1)*bsz2; j++)
                magOffDiagBlk += abs(sys.Mg[bsz2*pj_u+j]);
            for (Int j=0; j<nrj_l*bsz2; j++)
                magOffDiagBlk += abs(sys.Mg[bsz2*pj_l+j]);
            
            // Compute condition number (in infinity norm) of the diagonal block:
            double infNormHii = 0.0, infNormHiiInv = 0.0;
            for (Int i=0; i<bsz; i++) {
                double rowSum = 0.0;
                for (Int j=0; j<bsz; j++)
                    rowSum += abs(sys.Mg[bsz2*pj+i*bsz+j]);
                if (rowSum > infNormHii)
                    infNormHii = rowSum;
                
                double rowSumInv = 0.0;
                for (Int j=0; j<bsz; j++)
                    rowSumInv += abs(sys.HiiInvs[this_thread][i*bsz+j]);
                if (rowSumInv > infNormHiiInv)
                    infNormHiiInv = rowSumInv;
            }
            double inftyCondNumber = infNormHii * infNormHiiInv;
            
            if (inftyCondNumber > minInftyCondNumber) {
                for (Int j=0; j<bsz; j++) {
                    double entry = sys.Mg[bsz2*pj+j*(bsz+1)];
                    if (entry != 0.0)
                        sys.Mg[bsz2*pj+j*(bsz+1)] += (entry / abs(entry)) * infNormHii / (double) bsz;
                    else
                        sys.Mg[bsz2*pj+j*(bsz+1)] += infNormHii / (double) bsz;
                }
            }
            
//             if (inftyCondNumber > 1.0e4)
//                 printf("Condition number %g\n", inftyCondNumber);
            
// // //             if ((magDiagBlk/magOffDiagBlk) < epsDiagDom) {
// // //                 double multFactor = epsDiagDom / (magDiagBlk/magOffDiagBlk);
// // //                 for (Int j=0; j<bsz2; j++)
// // //                     sys.Mg[bsz2*pj+j] *= multFactor;
// // //                 
// // // //                 for (Int j=0; j<bsz2; j++) {
// // // //                     double entry = sys.Mg[bsz2*pj+j];
// // // //                     if (entry != 0) {
// // // //                         double sign_entry = entry / abs(entry);
// // // //                         sys.Mg[bsz2*pj+j] = epsDiagDom * sign_entry * magOffDiagBlk;
// // // //                     }
// // // //                 }
// // //             }
        }
    }
    
    // Compute BILUk preconditioner
    for (Int rj=0; rj<nrowsP; rj++) {
        Int rjj = sys.ordered2unordered[rj];
        
        Int pj = sys.Uoent2entStart[rj]-1;
        Int pj_l = sys.Loent2entStart[rj];
        Int pj_u = sys.Uoent2entStart[rj+1];
        
        // Compute inverse of block diagonal of Hg using LU decomposition
        DGETRF(&bsz,&bsz,&sys.Mg[bsz2*pj],&bsz,ipiv,&info);
        DGETRI(&bsz,&sys.Mg[bsz2*pj],&bsz,ipiv,work,&lwork,&info);
        
        Int nrjP = rowptsP[rjj+1]-rowptsP[rjj];
        Int nrj_l = sys.Loent2entStart[rj+1]-sys.Loent2entStart[rj];        // Neighbors in L
        Int nrj_u = sys.Uoent2entStart[rj]-sys.Uoent2entStart[rj+1];        // Neighbors in U (including the diagonal block)
        for (Int i=0; i<(nrj_u-1); i++) {
            Int pji_u = pj_u+i;
            Int rii = sys.LUoent2ent[pji_u];
            Int ri = sys.unordered2ordered[rii];
            
            Int pi = sys.Uoent2entStart[ri]-1;
            Int pi_l = sys.Loent2entStart[ri];
            Int pi_u = sys.Uoent2entStart[ri+1];
            
            Int nriP = rowptsP[rii+1]-rowptsP[rii];
            Int nri_l = sys.Loent2entStart[ri+1]-sys.Loent2entStart[ri];
            Int nri_u = sys.Uoent2entStart[ri]-sys.Uoent2entStart[ri+1];
            Int l;
            for (l=0; l<nri_l; l++) {
                if (sys.LUoent2ent[pi_l+l] == rjj) break;
            }
            if (l==nri_l && l!=0) {
                error("Error No. 3 in computeBILUk\n");
            }
            Int pij_l = pi_l + l;
            
            // Lij <- Uij*inv(Ujj)
            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.Mg[bsz2*pij_l], &bsz, &sys.Mg[bsz2*pj],
                  &bsz, &zero, Lij, &bsz);
            DCOPY(&bsz2, Lij, &inc, &sys.Mg[bsz2*pij_l], &inc);
            
            // Uii <- Uii - Lij*Uji
            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pij_l], &bsz, &sys.Mg[bsz2*pji_u],
                  &bsz, &one, &sys.Mg[bsz2*pi], &bsz);
            
            // L_il <- L_il - L_ij*U_jm (m "=" l)
            for (l=0; l<nri_l; l++) {
                Int pil_l = pi_l+l;
                Int rll  = sys.LUoent2ent[pil_l];
                Int rl = sys.unordered2ordered[rll];
                
                if (rl > rj) {
                    for (Int m=0; m<(nrj_u-1); m++) {
                        if (rll == sys.LUoent2ent[pj_u+m]) {
                            Int pjm_u = pj_u+m;
                            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pij_l], &bsz, &sys.Mg[bsz2*pjm_u],
                                  &bsz, &one, &sys.Mg[bsz2*pil_l], &bsz);
                        }
                    }
                }
            }
            
            // U_il <- U_il - L_ij*U_jm (m "=" l)
            for (l=0; l<(nri_u-1); l++) {
                Int pil_u = pi_u+l;
                Int rll  = sys.LUoent2ent[pil_u];
                Int rl = sys.unordered2ordered[rll];

                if (rl > rj) {
                    for (Int m=0; m<(nrj_u-1); m++) {
                        if (rll == sys.LUoent2ent[pj_u+m]) {
                            Int pjm_u = pj_u+m;
                            DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pij_l], &bsz, &sys.Mg[bsz2*pjm_u],
                                  &bsz, &one, &sys.Mg[bsz2*pil_u], &bsz);
                        }
                    }
                }
            }
        }
    }
}

// // // void computeBILU0_v2(sysstruct &sys, Int preconditioner)
// // // {
// // //     error("computeBILU0_v2 function is deprecated.\n");
// // //     // The block format of Mg matches the block format of Hg (rows and columns are stored in non-ordered way)
// // // 
// // //     Int nrows, nblks;
// // //     Int *rowpts, *colind;
// // //     if (preconditioner == 0) {
// // //         nrows = sys.numEntities;
// // //         nblks = sys.numBlocks;
// // //         rowpts = &sys.ent2entStart[0];
// // //         colind = &sys.ent2ent[0];
// // //     }
// // //     else {
// // //         nrows = sys.BJ_nrows;
// // //         nblks = sys.BJ_nblks;
// // //         rowpts = &sys.BJ_rowpts[0];
// // //         colind = &sys.BJ_colind[0];
// // //     }
// // // 
// // //     Int bsz = sys.blkSize;
// // //     Int bsz2 = bsz*bsz;
// // //     Int inc = 1, i,j,k,l,m,n, rj, jj, ri, pji, nri, nrj, pil, nnz;
// // //     Int info, ij, ii, ji, rl, il, jm;
// // // 
// // //     char chn = 'N';
// // //     double one = 1.0, zero = 0.0, minusone = -1.0;
// // // 
// // //     Int *ipiv = &sys.ipiv[0];
// // //     Int lwork = bsz;
// // //     double *work = &sys.work[0];
// // //     double *Lij = &sys.HiiInv[0];
// // // 
// // //     // Copy Hg into Mg with the desired storage format
// // //     nnz = sys.Mg.size();
// // //     DCOPY(&nnz, &sys.Hg[0], &inc, &sys.Mg[0], &inc);
// // // 
// // //     // Compute BILU0 preconditioner
// // //     for (j=0; j<nrows; j++) {
// // //         rj = sys.ordered2unordered[j];  // row rj
// // //         jj = rowpts[rj];
// // //         
// // //         // Compute inverse of block diagonal of Hg using LU decomposition
// // //         DGETRF(&bsz,&bsz,&sys.Mg[bsz2*jj],&bsz,ipiv,&info);
// // //         DGETRI(&bsz,&sys.Mg[bsz2*jj],&bsz,ipiv,work,&lwork,&info);
// // // 
// // //         nrj = rowpts[rj+1]-rowpts[rj];
// // //         for (i=1; i<nrj; i++) {
// // //             pji = jj+i;
// // //             ri = colind[pji];
// // //             if (sys.unordered2ordered[ri] > j) {
// // //                 ii = rowpts[ri];
// // //                 nri = rowpts[ri+1]-rowpts[ri];
// // //                 for (l = 0; l < nri; l++) {
// // //                     pil = ii+l;
// // //                     if (colind[pil] == rj) break;
// // //                 }
// // //                 if (l==nri)
// // //                     error("Error No. 1 in computeBILU0_v2\n");
// // //                 ij = pil;
// // // 
// // //                 // Lij <- Uij*inv(Ujj)
// // //                 DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &one, &sys.Mg[bsz2*ij], &bsz, &sys.Mg[bsz2*jj],
// // //                       &bsz, &zero, Lij, &bsz);
// // //                 DCOPY(&bsz2, Lij, &inc, &sys.Mg[bsz2*ij], &inc);
// // // 
// // //                 // Uii <- Uii - Lij*Uji
// // //                 ji = jj+i;
// // //                 DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &minusone, &sys.Mg[bsz2*ij], &bsz, &sys.Mg[bsz2*ji],
// // //                       &bsz, &one, &sys.Mg[bsz2*ii], &bsz);
// // // 
// // //                 for (l=1; l<nri; l++) {
// // //                     pil = ii+l;
// // //                     rl  = colind[pil];
// // //                     if (sys.unordered2ordered[rl] > j) {
// // //                         for (m=1; m<nrj; m++)
// // //                             if (rl == colind[jj+m]) {
// // //                                 il = ii+l;
// // //                                 jm = jj+m;
// // //                                 // U_il <- U_il - L_ij*U_jm
// // //                                 DGEMM(&chn, &chn, &bsz, &bsz, &bsz, &minusone, &sys.Mg[bsz2*ij], &bsz, &sys.Mg[bsz2*jm],
// // //                                       &bsz, &one, &sys.Mg[bsz2*il], &bsz);
// // //                             }
// // //                     }
// // //                 }
// // //             }
// // //         }
// // //     }
// // // }

void computeEntityInvBJ(sysstruct &sys)
{
    Int nrowsP = sys.nrowsP, bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int inc = 1, info, pj;
    
    Int *rowptsJ;
    rowptsJ = &sys.ent2entStart[0];
    
    char chn = 'N';
    double one = 1.0, zero = 0.0, minusone = -1.0;
    
    vector<vector<Int> >* ipivs = &sys.ipivs;
    Int lwork = bsz;
    vector<vector<double> >* works = &sys.works;
    
//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for private(pj, info)
        for (Int rj=0; rj<nrowsP; rj++) {
            pj = rowptsJ[rj];
            
            // Copy diagonal block
            DCOPY(&bsz2, &sys.Hg[bsz2*pj], &inc, &sys.Mg[bsz2*rj], &inc);

            // Compute inverse of block diagonal of Hg using LU decomposition
            DGETRF(&bsz,&bsz,&sys.Mg[bsz2*rj],&bsz,&ipivs[0][this_thread][0],&info);
            DGETRI(&bsz,&sys.Mg[bsz2*rj],&bsz,&ipivs[0][this_thread][0],&works[0][this_thread][0],&lwork,&info);
        }
//     }
}

//void applyBILU0_v2(sysstruct &sys, double * r, Int preconditioner)
//{
//    // The block format of Mg matches the block format of Hg (rows and columns are stored in non-ordered way)
//
//    Int nrows, nblks;
//    Int *rowpts, *colind;
//    if (preconditioner == 0) {
//        nrows = sys.numEntities;
//        nblks = sys.numBlocks;
//        rowpts = &sys.ent2entStart[0];
//        colind = &sys.ent2ent[0];
//    }
//    else {
//        nrows = sys.BJ_nrows;
//        nblks = sys.BJ_nblks;
//        rowpts = &sys.BJ_rowpts[0];
//        colind = &sys.BJ_colind[0];
//    }
//
//    Int bsz = sys.blkSize;
//    Int bsz2 = bsz*bsz;
//
//    char chn = 'N';
//    double one = 1.0, zero = 0.0, minusone = -1.0;
//
//    Int inc = 1, oneInt = 1, i,j,k,l,m,n, rj, jj, ji, ri, pji, nrj, info;
//
//    // Forward solve r <- L \ r
//    for (j=0; j<nrows; j++) {
//        rj = sys.ordered2unordered[j];  // row rj
//        nrj = rowpts[rj+1]-rowpts[rj];
//        // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
//        for (i=1; i<nrj; i++) {
//            pji = rowpts[rj]+i;
//            ri = colind[pji];
//            if (sys.unordered2ordered[ri] < j) {
//                DGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pji], &bsz, &r[bsz*ri],
//                      &inc, &one, &r[bsz*rj], &inc);
//            }
//        }
//        /* r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks */
//    }
//
//    // Backward solve r <- U \ r
//    for (j=0; j<nrows; j++) {
//        rj = sys.ordered2unordered[nrows-j-1];  // row rj
//        jj = rowpts[rj];
//        nrj = rowpts[rj+1]-rowpts[rj];
//        // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
//        for (i=1; i<nrj; i++) {
//            pji = rowpts[rj]+i;
//            ri = colind[pji];
//            if (sys.unordered2ordered[ri] > sys.unordered2ordered[rj]) {
//                DGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pji], &bsz, &r[bsz*ri],
//                      &inc, &one, &r[bsz*rj], &inc);
//            }
//        }
//        // Compute r_j <- U_jj \ r_j
//        DGEMV(&chn, &bsz, &bsz, &one, &sys.Mg[bsz2*jj], &bsz, &r[bsz*rj],
//              &inc, &zero, &sys.rtmp[0], &inc);
//        DCOPY(&bsz, &sys.rtmp[0], &inc, &r[bsz*rj], &inc);
//    }
//}

void applyBILU0_sp(sysstruct &sys, float * r_sp, double * r, Int preconditioner, Int * requestCounter)
{
    // Mg is stored in memory in the same order as is required for the BILU0 solve (L first, then U).
    
#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && preconditioner == 0 && (sys.reorderMethod == 0 || sys.reorderMethod == 1 || sys.reorderMethod == 2)) {
        recvvector(sys, &r[0], requestCounter);
    }
#endif
    
    if (sys.reorderMethod == 3 || sys.reorderMethod == 4) {
        error("Mixed-precision algorithm not implemented for MDF with constraints.\n");
    }
    
    Int nrows, nrowsBJ, nblks;
    Int *rowpts, *colind;
    float *rDenseRow_sp = &sys.rDenseRow_sp[0];
    
//     if (preconditioner == 0) {
//         nrows = sys.numEntities;
//         nblks = sys.numBlocks;
//         rowpts = &sys.ent2entStart[0];
//         colind = &sys.ent2ent[0];
//         if (sys.nproc > 1)
//             nrowsBJ = sys.BJ_nrows;
//         else
//             nrowsBJ = nrows;
//     }
    if (preconditioner == 1) {
        nrows = sys.BJ_nrows;
        nrowsBJ = nrows;
        nblks = sys.numBlocks;
        rowpts = &sys.ent2entStart[0];
        colind = &sys.ent2ent[0];
//         nrows = sys.BJ_nrows;
//         nrowsBJ = nrows;
//         nblks = sys.BJ_nblks;
//         rowpts = &sys.BJ_rowpts[0];
//         colind = &sys.BJ_colind[0];
    }
    else
        error("Preconditioner flag not recognized in applyBILU0_sp.\n");
    
    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int nb;

    char chn = 'N';
    float one = 1.0, zero = 0.0, minusone = -1.0;

    Int *Uoent2entStart = &sys.Uoent2entStart[0];
    Int *ordered2unordered = &sys.ordered2unordered[0];
    Int *LUoent2ent = &sys.LUoent2ent[0];
    Int *Loent2entStart = &sys.Loent2entStart[0];

    Int inc = 1, i,j,k,l,m,n,jjj;
    Int rj, rjj, ri, rii, pj, pj_l, pj_u, pjj, pjji, pji, pji_l, pji_u, nrj, nrj_l, nrj_u, info;

    for (i = 0; i < nrows*bsz; i++)
        r_sp[i] = (float) r[i];

    // Forward solve r <- L \ r
    if (sys.precSolveImplementation == 0) {
        for (rj=0; rj<nrowsBJ; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                SGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg_sp[bsz2*pji_l], &bsz, &r_sp[bsz*rii],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (rj=0; rj<nrowsBJ; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                for (jjj = 0; jjj < bsz; jjj++)
                   rDenseRow_sp[bsz*i+jjj] = r_sp[bsz*rii+jjj];
                //SCOPY(&bsz, &r_sp[bsz*rii], &inc, &rDenseRow_sp[bsz*i], &inc);
            }

            if (nrj_l > 0) {
                nb = bsz*nrj_l;
                SGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg_sp[bsz2*pj_l], &bsz, &rDenseRow_sp[0],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    
#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && preconditioner == 0 && (sys.reorderMethod == 3 || sys.reorderMethod == 4)) {
        recvvector(sys, &r[0], requestCounter);
    }
#endif
                        
    if (sys.precSolveImplementation == 0) {
        for (rj=nrowsBJ; rj<nrows; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                SGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg_sp[bsz2*pji_l], &bsz, &r_sp[bsz*rii],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (rj=nrowsBJ; rj<nrows; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                for (jjj = 0; jjj < bsz; jjj++)
                    rDenseRow_sp[bsz*i+jjj] = r_sp[bsz*rii+jjj];
//                 SCOPY(&bsz, &r_sp[bsz*rii], &inc, &rDenseRow_sp[bsz*i], &inc);
            }

            if (nrj_l > 0) {
                nb = bsz*nrj_l;
                SGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg_sp[bsz2*pj_l], &bsz, &rDenseRow_sp[0],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
                    
    // Backward solve r <- U \ r
    if (sys.precSolveImplementation == 0) {
        for (j=0; j<nrows; j++) {
            rj = nrows-j-1;
            rjj = ordered2unordered[rj];

            pj = Uoent2entStart[rj]-1;
            pj_u = Uoent2entStart[rj+1];

            nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];
 
            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            for (i=0; i<(nrj_u-1); i++) {
                pji_u = pj_u+i;
                rii = LUoent2ent[pji_u];
                SGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg_sp[bsz2*pji_u], &bsz, &r_sp[bsz*rii],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }

            // Compute r_j <- U_jj \ r_j
            SGEMV(&chn, &bsz, &bsz, &one, &sys.Mg_sp[bsz2*pj], &bsz, &r_sp[bsz*rjj],
                  &inc, &zero, &sys.rtmp_sp[0], &inc);
            
            SCOPY(&bsz, &sys.rtmp_sp[0], &inc, &r_sp[bsz*rjj], &inc);
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (j=0; j<nrows; j++) {
            rj = nrows-j-1;
            rjj = ordered2unordered[rj];

            pj = Uoent2entStart[rj]-1;
            pj_u = Uoent2entStart[rj+1];

            nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];

            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            for (i=0; i<(nrj_u-1); i++) {
                pji_u = pj_u+i;
                rii = LUoent2ent[pji_u];
                for (jjj = 0; jjj < bsz; jjj++)
                    rDenseRow_sp[bsz*i+jjj] = r_sp[bsz*rii+jjj];
                //SCOPY(&bsz, &r_sp[bsz*rii], &inc, &rDenseRow_sp[bsz*i], &inc);
            }

            if (nrj_u > 1) {
                nb = bsz*(nrj_u-1);
                SGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg_sp[bsz2*pj_u], &bsz, &rDenseRow_sp[0],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }

            // Compute r_j <- U_jj \ r_j
            SGEMV(&chn, &bsz, &bsz, &one, &sys.Mg_sp[bsz2*pj], &bsz, &r_sp[bsz*rjj],
                  &inc, &zero, &sys.rtmp_sp[0], &inc);
            SCOPY(&bsz, &sys.rtmp_sp[0], &inc, &r_sp[bsz*rjj], &inc);
        }
    }
                    
    for (Int i = 0; i < nrowsBJ*bsz; i++)
        r[i] = (double) r_sp[i];
}

void applyBILU0_dp(sysstruct &sys, double * r, Int preconditioner, Int * requestCounter)
{
    // Mg is stored in memory in the same order as is required for the BILU0 solve (L first, then U).
    
#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && preconditioner == 0 && (sys.reorderMethod == 0 || sys.reorderMethod == 1 || sys.reorderMethod == 2)) {
        recvvector(sys, &r[0], requestCounter);
    }
#endif
    
    Int nrows, nrowsBJ, nblks;
    Int *rowpts, *colind;
    double *rDenseRow = &sys.xDenseRow[0][0];
    
//     if (preconditioner == 0) {
//         nrows = sys.numEntities;
//         nblks = sys.numBlocks;
//         rowpts = &sys.ent2entStart[0];
//         colind = &sys.ent2ent[0];
//         if (sys.nproc > 1)
//             nrowsBJ = sys.BJ_nrows;
//         else
//             nrowsBJ = nrows;
//     }
    if (preconditioner == 1) {
        nrows = sys.BJ_nrows;
        nrowsBJ = nrows;
        nblks = sys.numBlocks;
        rowpts = &sys.ent2entStart[0];
        colind = &sys.ent2ent[0];
//         nrows = sys.BJ_nrows;
//         nrowsBJ = nrows;
//         nblks = sys.BJ_nblks;
//         rowpts = &sys.BJ_rowpts[0];
//         colind = &sys.BJ_colind[0];
    }
    else
        error("Preconditioner flag not recognized in applyBILU0_dp.\n");

    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int nb;

    char chn = 'N';
    double one = 1.0, zero = 0.0, minusone = -1.0;

    Int *Uoent2entStart = &sys.Uoent2entStart[0];
    Int *ordered2unordered = &sys.ordered2unordered[0];
    Int *LUoent2ent = &sys.LUoent2ent[0];
    Int *Loent2entStart = &sys.Loent2entStart[0];

    Int inc = 1, i,j,k,l,m,n,jjj;
    Int rj, rjj, ri, rii, pj, pj_l, pj_u, pjj, pjji, pji, pji_l, pji_u, nrj, nrj_l, nrj_u, info;

    // Forward solve r <- L \ r
    if (sys.precSolveImplementation == 0) {
        for (rj=0; rj<nrowsBJ; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                DGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pji_l], &bsz, &r[bsz*rii],
                      &inc, &one, &r[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (rj=0; rj<nrowsBJ; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                for (jjj = 0; jjj < bsz; jjj++)
                   rDenseRow[bsz*i+jjj] = r[bsz*rii+jjj];
                //DCOPY(&bsz, &r[bsz*rii], &inc, &rDenseRow[bsz*i], &inc);
            }

            if (nrj_l > 0) {
                nb = bsz*nrj_l;
                DGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg[bsz2*pj_l], &bsz, &rDenseRow[0],
                      &inc, &one, &r[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }

#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && preconditioner == 0 && (sys.reorderMethod == 3 || sys.reorderMethod == 4)) {
        recvvector(sys, &r[0], requestCounter);
    }
#endif

    if (sys.precSolveImplementation == 0) {
        for (rj=nrowsBJ; rj<nrows; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                DGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pji_l], &bsz, &r[bsz*rii],
                      &inc, &one, &r[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (rj=nrowsBJ; rj<nrows; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                for (jjj = 0; jjj < bsz; jjj++)
                   rDenseRow[bsz*i+jjj] = r[bsz*rii+jjj];
                //DCOPY(&bsz, &r[bsz*rii], &inc, &rDenseRow[bsz*i], &inc);
            }

            if (nrj_l > 0) {
                nb = bsz*nrj_l;
                DGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg[bsz2*pj_l], &bsz, &rDenseRow[0],
                      &inc, &one, &r[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }

    // Backward solve r <- U \ r
    if (sys.precSolveImplementation == 0) {
        for (j=0; j<nrows; j++) {
            rj = nrows-j-1;
            rjj = ordered2unordered[rj];

            pj = Uoent2entStart[rj]-1;
            pj_u = Uoent2entStart[rj+1];

            nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];

            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            for (i=0; i<(nrj_u-1); i++) {
                pji_u = pj_u+i;
                rii = LUoent2ent[pji_u];
                DGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pji_u], &bsz, &r[bsz*rii],
                      &inc, &one, &r[bsz*rjj], &inc);
            }
            
            // Compute r_j <- U_jj \ r_j
            DGEMV(&chn, &bsz, &bsz, &one, &sys.Mg[bsz2*pj], &bsz, &r[bsz*rjj],
                  &inc, &zero, &sys.rtmp[0], &inc);
            DCOPY(&bsz, &sys.rtmp[0], &inc, &r[bsz*rjj], &inc);
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (j=0; j<nrows; j++) {
            rj = nrows-j-1;
            rjj = ordered2unordered[rj];

            pj = Uoent2entStart[rj]-1;
            pj_u = Uoent2entStart[rj+1];

            nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];

            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            for (i=0; i<(nrj_u-1); i++) {
                pji_u = pj_u+i;
                rii = LUoent2ent[pji_u];
                for (jjj = 0; jjj < bsz; jjj++)
                   rDenseRow[bsz*i+jjj] = r[bsz*rii+jjj];
                //DCOPY(&bsz, &r[bsz*rii], &inc, &rDenseRow[bsz*i], &inc);
            }

            if (nrj_u > 1) {
                nb = bsz*(nrj_u-1);
                DGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg[bsz2*pj_u], &bsz, &rDenseRow[0],
                      &inc, &one, &r[bsz*rjj], &inc);
            }

            // Compute r_j <- U_jj \ r_j
            DGEMV(&chn, &bsz, &bsz, &one, &sys.Mg[bsz2*pj], &bsz, &r[bsz*rjj],
                  &inc, &zero, &sys.rtmp[0], &inc);
            DCOPY(&bsz, &sys.rtmp[0], &inc, &r[bsz*rjj], &inc);
        }
    }
}

void applyBILUk_sp(sysstruct &sys, float * r_sp, double * r, Int preconditioner, Int * requestCounter)
{
    // Mg is stored in memory in the same order as is required for the BILUk solve (L first, then U).
    
#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && preconditioner == 0 && (sys.reorderMethod == 0 || sys.reorderMethod == 1 || sys.reorderMethod == 2)) {
        recvvector(sys, &r[0], requestCounter);
    }
#endif
    
    if (sys.reorderMethod == 3 || sys.reorderMethod == 4) {
        error("Mixed-precision algorithm not implemented for MDF with constraints.\n");
    }
    
    char chn = 'N';
    float one = 1.0, zero = 0.0, minusone = -1.0;
    
    float *rDenseRow_sp = &sys.rDenseRow_sp[0];
    
    Int nrowsP, nrowsBJ;
    Int nblks = sys.numBlocks;
    if (preconditioner == 0 || preconditioner == 1) {
        nrowsP = sys.nrowsP;
        if (sys.nproc > 1)
            nrowsBJ = sys.BJ_nrows;
        else
            nrowsBJ = nrowsP;
    }
    else
        error("Preconditioner flag not recognized in applyBILUk_sp.\n");
    
    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int nb;
    
    Int *Uoent2entStart = &sys.Uoent2entStart[0];
    Int *ordered2unordered = &sys.ordered2unordered[0];
    Int *LUoent2ent = &sys.LUoent2ent[0];
    Int *Loent2entStart = &sys.Loent2entStart[0];
    
    Int inc = 1, i,j,k,l,m,n,jjj;
    Int rj, rjj, ri, rii, pj, pj_l, pj_u, pjj, pjji, pji, pji_l, pji_u, nrj, nrj_l, nrj_u, info;

    for (i = 0; i < nrowsP*bsz; i++)
        r_sp[i] = (float) r[i];

    // Forward solve r <- L \ r
    if (sys.precSolveImplementation == 0) {
        for (rj=0; rj<nrowsBJ; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                SGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg_sp[bsz2*pji_l], &bsz, &r_sp[bsz*rii],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (rj=0; rj<nrowsBJ; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                for (jjj = 0; jjj < bsz; jjj++)
                   rDenseRow_sp[bsz*i+jjj] = r_sp[bsz*rii+jjj];
                //SCOPY(&bsz, &r_sp[bsz*rii], &inc, &rDenseRow_sp[bsz*i], &inc);
            }

            if (nrj_l > 0) {
                nb = bsz*nrj_l;
                SGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg_sp[bsz2*pj_l], &bsz, &rDenseRow_sp[0],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    
#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && preconditioner == 0 && (sys.reorderMethod == 3 || sys.reorderMethod == 4)) {
        recvvector(sys, &r[0], requestCounter);
    }
#endif
                        
    if (sys.precSolveImplementation == 0) {
        for (rj=nrowsBJ; rj<nrowsP; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                SGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg_sp[bsz2*pji_l], &bsz, &r_sp[bsz*rii],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (rj=nrowsBJ; rj<nrowsP; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                for (jjj = 0; jjj < bsz; jjj++)
                    rDenseRow_sp[bsz*i+jjj] = r_sp[bsz*rii+jjj];
//                 SCOPY(&bsz, &r_sp[bsz*rii], &inc, &rDenseRow_sp[bsz*i], &inc);
            }

            if (nrj_l > 0) {
                nb = bsz*nrj_l;
                SGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg_sp[bsz2*pj_l], &bsz, &rDenseRow_sp[0],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
                    
    // Backward solve r <- U \ r
    if (sys.precSolveImplementation == 0) {
        for (j=0; j<nrowsP; j++) {
            rj = nrowsP-j-1;
            rjj = ordered2unordered[rj];

            pj = Uoent2entStart[rj]-1;
            pj_u = Uoent2entStart[rj+1];

            nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];
 
            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            for (i=0; i<(nrj_u-1); i++) {
                pji_u = pj_u+i;
                rii = LUoent2ent[pji_u];
                SGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg_sp[bsz2*pji_u], &bsz, &r_sp[bsz*rii],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }

            // Compute r_j <- U_jj \ r_j
            SGEMV(&chn, &bsz, &bsz, &one, &sys.Mg_sp[bsz2*pj], &bsz, &r_sp[bsz*rjj],
                  &inc, &zero, &sys.rtmp_sp[0], &inc);
            
            SCOPY(&bsz, &sys.rtmp_sp[0], &inc, &r_sp[bsz*rjj], &inc);
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (j=0; j<nrowsP; j++) {
            rj = nrowsP-j-1;
            rjj = ordered2unordered[rj];

            pj = Uoent2entStart[rj]-1;
            pj_u = Uoent2entStart[rj+1];

            nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];

            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            for (i=0; i<(nrj_u-1); i++) {
                pji_u = pj_u+i;
                rii = LUoent2ent[pji_u];
                for (jjj = 0; jjj < bsz; jjj++)
                    rDenseRow_sp[bsz*i+jjj] = r_sp[bsz*rii+jjj];
                //SCOPY(&bsz, &r_sp[bsz*rii], &inc, &rDenseRow_sp[bsz*i], &inc);
            }

            if (nrj_u > 1) {
                nb = bsz*(nrj_u-1);
                SGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg_sp[bsz2*pj_u], &bsz, &rDenseRow_sp[0],
                      &inc, &one, &r_sp[bsz*rjj], &inc);
            }

            // Compute r_j <- U_jj \ r_j
            SGEMV(&chn, &bsz, &bsz, &one, &sys.Mg_sp[bsz2*pj], &bsz, &r_sp[bsz*rjj],
                  &inc, &zero, &sys.rtmp_sp[0], &inc);
            SCOPY(&bsz, &sys.rtmp_sp[0], &inc, &r_sp[bsz*rjj], &inc);
        }
    }
                    
    for (Int i = 0; i < nrowsBJ*bsz; i++)
        r[i] = (double) r_sp[i];
}

void applyBILUk_dp(sysstruct &sys, double * r, Int preconditioner, Int * requestCounter)
{
    // Mg is stored in memory in the same order as is required for the BILUk solve (L first, then U).

#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && preconditioner == 0 && (sys.reorderMethod == 0 || sys.reorderMethod == 1 || sys.reorderMethod == 2)) {
        recvvector(sys, &r[0], requestCounter);
    }
#endif
    
    char chn = 'N';
    double one = 1.0, zero = 0.0, minusone = -1.0;
    
    double *rDenseRow = &sys.xDenseRow[0][0];
    
    Int nrowsP, nrowsBJ;
    Int nblks = sys.numBlocks;
    
    if (preconditioner == 0 || preconditioner == 1) {
        nrowsP = sys.nrowsP;
        if (sys.nproc > 1)
            nrowsBJ = sys.BJ_nrows;
        else
            nrowsBJ = nrowsP;
    }
    else
        error("Preconditioner flag not recognized in applyBILUk_sp.\n");
    
    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int nb;
    
    Int *Uoent2entStart = &sys.Uoent2entStart[0];
    Int *ordered2unordered = &sys.ordered2unordered[0];
    Int *LUoent2ent = &sys.LUoent2ent[0];
    Int *Loent2entStart = &sys.Loent2entStart[0];

    Int inc = 1, i,j,k,l,m,n,jjj;
    Int rj, rjj, ri, rii, pj, pj_l, pj_u, pjj, pjji, pji, pji_l, pji_u, nrj, nrj_l, nrj_u, info;

    // Forward solve r <- L \ r
    if (sys.precSolveImplementation == 0) {
        for (rj=0; rj<nrowsBJ; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                DGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pji_l], &bsz, &r[bsz*rii],
                      &inc, &one, &r[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (rj=0; rj<nrowsBJ; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                for (jjj = 0; jjj < bsz; jjj++)
                   rDenseRow[bsz*i+jjj] = r[bsz*rii+jjj];
                //DCOPY(&bsz, &r[bsz*rii], &inc, &rDenseRow[bsz*i], &inc);
            }

            if (nrj_l > 0) {
                nb = bsz*nrj_l;
                DGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg[bsz2*pj_l], &bsz, &rDenseRow[0],
                      &inc, &one, &r[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }

#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && preconditioner == 0 && (sys.reorderMethod == 3 || sys.reorderMethod == 4)) {
        recvvector(sys, &r[0], requestCounter);
    }
#endif

    if (sys.precSolveImplementation == 0) {
        for (rj=nrowsBJ; rj<nrowsP; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                DGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pji_l], &bsz, &r[bsz*rii],
                      &inc, &one, &r[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (rj=nrowsBJ; rj<nrowsP; rj++) {
            rjj = ordered2unordered[rj];
            pj_l = Loent2entStart[rj];
            nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            for (i=0; i<nrj_l; i++) {
                pji_l = pj_l+i;
                rii = LUoent2ent[pji_l];
                for (jjj = 0; jjj < bsz; jjj++)
                   rDenseRow[bsz*i+jjj] = r[bsz*rii+jjj];
                //DCOPY(&bsz, &r[bsz*rii], &inc, &rDenseRow[bsz*i], &inc);
            }

            if (nrj_l > 0) {
                nb = bsz*nrj_l;
                DGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg[bsz2*pj_l], &bsz, &rDenseRow[0],
                      &inc, &one, &r[bsz*rjj], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }

    // Backward solve r <- U \ r
    if (sys.precSolveImplementation == 0) {
        for (j=0; j<nrowsP; j++) {
            rj = nrowsP-j-1;
            rjj = ordered2unordered[rj];

            pj = Uoent2entStart[rj]-1;
            pj_u = Uoent2entStart[rj+1];

            nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];

            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            for (i=0; i<(nrj_u-1); i++) {
                pji_u = pj_u+i;
                rii = LUoent2ent[pji_u];
                DGEMV(&chn, &bsz, &bsz, &minusone, &sys.Mg[bsz2*pji_u], &bsz, &r[bsz*rii],
                      &inc, &one, &r[bsz*rjj], &inc);
            }

            // Compute r_j <- U_jj \ r_j
            DGEMV(&chn, &bsz, &bsz, &one, &sys.Mg[bsz2*pj], &bsz, &r[bsz*rjj],
                  &inc, &zero, &sys.rtmp[0], &inc);
            DCOPY(&bsz, &sys.rtmp[0], &inc, &r[bsz*rjj], &inc);
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (j=0; j<nrowsP; j++) {
            rj = nrowsP-j-1;
            rjj = ordered2unordered[rj];

            pj = Uoent2entStart[rj]-1;
            pj_u = Uoent2entStart[rj+1];

            nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];

            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            for (i=0; i<(nrj_u-1); i++) {
                pji_u = pj_u+i;
                rii = LUoent2ent[pji_u];
                for (jjj = 0; jjj < bsz; jjj++)
                   rDenseRow[bsz*i+jjj] = r[bsz*rii+jjj];
                //DCOPY(&bsz, &r[bsz*rii], &inc, &rDenseRow[bsz*i], &inc);
            }

            if (nrj_u > 1) {
                nb = bsz*(nrj_u-1);
                DGEMV(&chn, &bsz, &nb, &minusone, &sys.Mg[bsz2*pj_u], &bsz, &rDenseRow[0],
                      &inc, &one, &r[bsz*rjj], &inc);
            }

            // Compute r_j <- U_jj \ r_j
            DGEMV(&chn, &bsz, &bsz, &one, &sys.Mg[bsz2*pj], &bsz, &r[bsz*rjj],
                  &inc, &zero, &sys.rtmp[0], &inc);
            DCOPY(&bsz, &sys.rtmp[0], &inc, &r[bsz*rjj], &inc);
        }
    }
}


void applyEntityInvBJ_sp(sysstruct &sys, float * r_sp, double * r)
{
    Int nrowsP = sys.nrowsP, bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int inc = 1, pj;
    
    char chn = 'N';
    float one = 1.0, zero = 0.0, minusone = -1.0;
    
//     #pragma omp parallel for num_threads(sys.noThreads)
    for (Int i = 0; i < nrowsP*bsz; i++)
        r_sp[i] = (float) r[i];
    
//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for private(pj)
        for (Int rj=0; rj<nrowsP; rj++) {
            SGEMV(&chn, &bsz, &bsz, &one, &sys.Mg_sp[bsz2*rj], &bsz, &r_sp[bsz*rj],
                  &inc, &zero, &sys.rtmps_sp[this_thread][0], &inc);
            SCOPY(&bsz, &sys.rtmps_sp[this_thread][0], &inc, &r_sp[bsz*rj], &inc);
        }
//     }
    
//     #pragma omp parallel for num_threads(sys.noThreads)
    for (Int i = 0; i < nrowsP*bsz; i++)
        r[i] = (double) r_sp[i];
}

void applyEntityInvBJ_dp(sysstruct &sys, double * r)
{
    Int nrowsP = sys.nrowsP, bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    Int inc = 1, pj;
    
    char chn = 'N';
    double one = 1.0, zero = 0.0, minusone = -1.0;
    
//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for private(pj)
        for (Int rj=0; rj<nrowsP; rj++) {
            DGEMV(&chn, &bsz, &bsz, &one, &sys.Mg[bsz2*rj], &bsz, &r[bsz*rj],
                  &inc, &zero, &sys.rtmps[this_thread][0], &inc);
            DCOPY(&bsz, &sys.rtmps[this_thread][0], &inc, &r[bsz*rj], &inc);
        }
//     }
}

void computeOrdering(sysstruct &sys, Int preconditioner)
{
    // TODO: Change the code so that for preconditioner == -1 and preconditioner == 2 no reorder is even computed
    if ((preconditioner == -1 || preconditioner == 2) && sys.reorderMethod != 0)
        sys.reorderMethod = 0;
    
    if (sys.reorderMethod == 0)             // No reorder
        computeNoOrdering(sys, preconditioner);
    else if (sys.reorderMethod == 1)        // Approximate MDF
        computeApproximateMDFordering(sys, preconditioner);
    else if (sys.reorderMethod == 2)        // Exact MDF
        computeExactMDFordering(sys, preconditioner);
    else if (sys.reorderMethod == 3)        // Approximate MDF with restricted order for entities not in subdomain (only for parallel RAS preconditioner)
        computeApproximateMDForderingWithConstraints(sys, preconditioner);
    else if (sys.reorderMethod == 4)        // Exact MDF with restricted order for entities not in subdomain (only for parallel RAS preconditioner)
        computeExactMDForderingWithConstraints(sys, preconditioner);
    else
        printf("Ordering algorithm not implemented yet.\n");
    
    if (preconditioner == 0)
        postprocessOrderingForBILUk(sys);
    else if (preconditioner == 1)
        postprocessOrderingForBILUk(sys);
}

void computePreconditioner(sysstruct &sys, Int preconditioner)
{
    if (preconditioner == -1) {
    }
    else if (preconditioner == 0)
        computeBILUk(sys);
    else if (preconditioner == 1)
        computeBILUk(sys);
    else if (preconditioner == 2)
        computeEntityInvBJ(sys);
    else
        error("Preconditioner flag not recognized in computePreconditioner.\n");
    
    // Convert preconditioner into single precision:
    if (preconditioner != -1 && sys.precPrecision == 0) {
        for (Int i = 0; i < sys.Mg.size(); i++)
            sys.Mg_sp[i] = (float) sys.Mg[i];
    }
}

void applyPreconditioner_sp(sysstruct &sys, float * r_sp, double * r, Int preconditioner, Int * requestCounter)
{   
    if (sys.preconditioner == -1) {
    }
    else if (sys.preconditioner == 0)
        applyBILUk_sp(sys, r_sp, r, preconditioner, requestCounter);
    else if (sys.preconditioner == 1)
        applyBILUk_sp(sys, r_sp, r, preconditioner, requestCounter);
    else if (sys.preconditioner == 2)
        applyEntityInvBJ_sp(sys, r_sp, r);
    else
        error("Preconditioner flag not recognized in applyPreconditioner_sp.\n");
}

void applyPreconditioner_dp(sysstruct &sys, double * r, Int preconditioner, Int * requestCounter)
{   
    if (sys.preconditioner == -1) {
    }
    else if (sys.preconditioner == 0)
        applyBILUk_dp(sys, r, preconditioner, requestCounter);
    else if (sys.preconditioner == 1)
        applyBILUk_dp(sys, r, preconditioner, requestCounter);
    else if (sys.preconditioner == 2)
        applyEntityInvBJ_dp(sys, r);
    else
        error("Preconditioner flag not recognized in applyPreconditioner_dp.\n");
}

void applyPreconditioner(sysstruct &sys, double * r, Int preconditioner, Int * requestCounter)
{   
    if (sys.precPrecision == 0 && sys.robustMode == 0) {
        float * r_sp = &sys.r_sp[0];
        applyPreconditioner_sp(sys, r_sp, r, preconditioner, requestCounter);
    }
    else
        applyPreconditioner_dp(sys, r, preconditioner, requestCounter);
}

#endif
