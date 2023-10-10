#ifndef __ORDERING
#define __ORDERING

// Written by: C. Nguyen & P. Fernandez

Int indexofSmallestElement(double* array, Int size)
{
    Int index = 0;
    
    for (Int i = 1; i < size; i++) {
        if (isnan(array[i]))
            error("NaN detected in indexofSmallestElement.\n");
        else if (isinf(array[i]))
            error("Inf detected in indexofSmallestElement.\n");
        else (array[i] < array[index])
            index = i;
    }
    return index;
}

// // // // void copyHg2MgWithHgPlusFillFormat(sysstruct &sys)
// // // // {
// // // //     Int pjiJ, inc = 1;
// // // //     Int bsz = sys.blkSize;
// // // //     Int bsz2 = bsz*bsz;
// // // //     
// // // //     Int nrowsP = sys.numEntitiesP;
// // // //     Int nblksP = sys.numBlocksP;
// // // //     Int *rowptsP = &sys.ent2entStartP[0];
// // // //     Int *rowptsJ = &sys.ent2entStartJ[0];
// // // //     Int *colindP = &sys.ent2entP[0];
// // // //     
// // // //     // Copy Hg into Mg with the desired storage format
// // // //     for (Int rj=0; rj<nrowsP; rj++) {
// // // //         Int pjJ = rowptsJ[rj];
// // // //         Int pjP = rowptsP[rj];
// // // //         
// // // //         Int nrjP = rowptsP[rj+1]-rowptsP[rj];
// // // //         
// // // //         // Copy diagonal block
// // // //         DCOPY(&bsz2, &sys.Hg[bsz2*pjJ], &inc, &sys.Mg[bsz2*pjP], &inc);
// // // //         
// // // //         // Copy non-diagonal blocks
// // // //         Int pjiJ = pjJ+1;
// // // //         for (Int i=1; i<nrjP; i++) {
// // // //             Int pjiP = pjP + i;
// // // //             
// // // //             if (! matchSparsity(rj,i)) {
// // // //                 for (j=0; j<sz2; j++)
// // // //                     sys.Mg[bsz2*pjiP+j] = 0.0;
// // // //             }
// // // //             else {
// // // //                 DCOPY(&bsz2, &sys.Hg[bsz2*pjiJ], &inc, &sys.Mg[bsz2*pjiP], &inc);
// // // //                 pjiJ++;
// // // //             }
// // // //         }
// // // //     }
// // // // }

void checkw(double *w, Int len, double bign)
{
    // Check sys.w vector:
    for (Int i = 0; i < len; i++) {
        if (isnan(w[i]) || isinf(w[i]) || (w[i] > bign)) {
            printf("Error in MDF reorder algorithm: w[%d] = %g.\nExecution will be terminated.\n", i, w[i]);
            exit(-1);
        }
    }
}

void computeNoOrdering(sysstruct &sys)
{
    for (Int i = 0; i < sys.numRowsP; i++) {
        sys.unordered2ordered[i] = i;
        sys.ordered2unordered[i] = i;
    }
}

void computeApproximateMDFordering(sysstruct &sys)
{
    /* TODO: Further optimization: If necessary, this function can be optimized by not recomputing w[i] when needs to be updated (we can just remove the contributions that no longer exist). */
    
    /* Note:
     * In this function, entity refers to faces in HDG and to face nodes in EDG and IEDG.
     * We assume all blocks are of the same size (i.e., nch (EDG) or nch*npf (UDH) are constant throughout the mesh).
     *      If not, an array indicating the size of each block and an array indicating the starting position each block in the matrix J would be required (easy fix).
     * The complexity of the MDF ordering is O(N^2) due to the indexofSmallestElement operation.
     * */
    
    Int inc = 1, info;
    Int nrowsP = sys.numEntitiesP;
    Int nblksP = sys.numBlocksP;
    Int *entSz = &sys.entSz[0];
    Int *ent2entStartP = &sys.ent2entStartP[0];
    Int *ent2entP = &sys.ent2entP[0];
    Int *rowStartP = &sys.rowStartP[0];
    double one = 1.0, zero = 0.0;
    double bign = 1.0e100;
    char chn = 'N';

    vector<vector<Int> >* ipivs = &sys.ipivs;
    vector<vector<double> >* works = &sys.works;
    Int *pl_local = new Int[sys.noThreads];
    
    //copyHg2MgWithHgPlusFillFormat(sys); ??
    copyHg2MgWithMgFormat(sys); ??
    
//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for //private(info, j, ii, nri, pij)
        for (Int ri=0; ri<nrowsP; ri++) {
            Int piSpP = ent2entStartP[ri];
            Int pijP = rowStartP[ri];
            Int nri = ent2entStartP[ri+1] - ent2entStartP[ri];
            Int bsz_row = entSz[ri];
            Int bsz2 = bsz_row * bsz_row;
            Int lwork = bsz_row;
            
            /* HiiInv <- inv(Hg_ii) */
            DCOPY(&bsz2, &sys.Mg[pijP], &inc, &sys.HiiInvs[this_thread][0], &inc);
            DGETRF(&bsz_row, &bsz_row, &sys.HiiInvs[this_thread][0], &bsz_row, &ipivs[0][this_thread][0], &info);
            DGETRI(&bsz_row, &sys.HiiInvs[this_thread][0], &bsz_row, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
            
            for (Int j=0; j<nri; j++) {
                Int pijSpP = piSpP+j;
                Int rj = ent2entP[pijSpP];
                Int bsz_col = entSz[rj];
                Int bsz2 = bsz_row * bsz_col;
                
                /* HiiInvHij <- HiiInv*Hg_ij (=inv(Hg_ii)*Hg_ij)  */
                if (matchSparsity(ri,j)) {
                    DGEMM(&chn, &chn, &bsz_row, &bsz_col, &bsz_row, &one, &sys.HiiInvs[this_thread][0], &bsz_row, &sys.Mg[pijP],
                          &bsz_col, &zero, &sys.HiiInvHijs[this_thread][0], &bsz_row);
                    sys.C[pijSpP] = DNRM2(&bsz2, &sys.HiiInvHijs[this_thread][0], &inc);
                }
                else {
                    sys.C[pijSpP] = 0.0;
                }
                pijP += bsz2;
            }
        }
        
//         #pragma omp for //private(i, j, m, ri, rj, nri, nrk, ik, kj, kk, pij, pim, pki, Cik, Ckj, shared)
        for (Int rk=0; rk<nrowsP; rk++) {
            sys.w[rk] = 0.0;
            Int pkSpP = ent2entStartP[rk];
            Int pkiP = rowStartP[rk];
            Int nrk = ent2entStartP[rk+1] - ent2entStartP[rk];
            for (Int i = 1; i<nrk; i++) {
                Int pkiSpP = pkSpP + i;
                Int ri = ent2entP[pkiSpP];
                Int nri = ent2entStartP[ri+1] - ent2entStartP[ri];
                Int pijSpP;
                for (Int j=1; j<nri; j++) {
                    pijSpP = ent2entStartP[ri]+j;
                    if (ent2entP[pijSpP] == rk) break;
                }
                if (j == nri)
                    error("Error G45BH0 in computeApproximateMDFordering.\n");
                Int pikSpP = pijSpP;
                double Cik = sys.C[pikSpP];
                for (Int j = 1; j<nrk; j++)
                    if (j != i) {
                        Int pkjSpP = pkSpP+j;
                        Int rj = ent2entP[pkjSpP];
                        double Ckj = sys.C[pkjSpP];

                        Int shared = 0;
                        for (Int m=1; m<nri; m++) {     // Determine if j neighbors i
                            Int pimSpP = ent2entStartP[ri]+m;
                            if (ent2entP[pimSpP] == rj) {
                                shared = 1;
                                break;
                            }
                        }
                        if (shared == 0)          // If j does not neighbor i, then increase the discarded fill in
                            sys.w[rk] += Cik*Cik*Ckj*Ckj;
                    }
            }
            sys.w[rk] = sqrt(sys.w[rk]);
        }
//     }
        
        // Check NaN and Inf in sys.w:
        checkw(&sys.w[0], nrowsP, bign);
        
        for (Int rl=0; rl<nrowsP; rl++) {
// // //             Int chunkLen, startPt;
// // //             splitArray1D(nrowsP, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //             pl_local[this_thread] = startPt + indexofSmallestElement(&sys.w[startPt], chunkLen);
            Int pl = indexofSmallestElement(&sys.w[0], nrowsP);
// // //             #pragma omp barrier     // TODO: This barrier (or in general the threaded way of computing pl may be too slow)
// // //             #pragma omp single
// // //             {
// // //                 Int pl = pl_local[0];
// // //                 for (Int i_th = 1; i_th < sys.noThreads; i_th++)
// // //                     if (sys.w[pl_local[i_th]] < sys.w[pl])
// // //                         pl = pl_local[i_th];
                sys.ordered2unordered[rl] = pl;
                sys.w[pl] = bign;
                Int nrl = ent2entStartP[pl+1] - ent2entStartP[pl];
// // //             }
// // //             #pragma omp for //private(i, j, k, m, ri, rj, rk, nri, nrk, pij, pki, plk, pim, ik, kj, Cik, Ckj, shared)
            for (Int k = 1; k<nrl; k++) {
                Int plkSpP = ent2entStartP[pl]+k;
                Int rk = ent2entP[plkSpP];
                if (sys.w[rk] < bign) {
                    sys.w[rk] = 0.0;
                    Int nrk = ent2entStartP[rk+1] - ent2entStartP[rk];
                    for (Int i = 1; i<nrk; i++) {
                        Int pkiSpP = ent2entStartP[rk]+i;
                        Int ri = ent2entP[pkiSpP];
                        if (sys.w[ri] < bign) {
                            Int nri = ent2entStartP[ri+1] - ent2entStartP[ri];
                            Int pijSpP;
                            for (Int j = 1; j<nri; j++) {
                                pijSpP = ent2entStartP[ri]+j;
                                if (ent2entP[pijSpP] == rk) break;
                            }
                            if (j == nri)
                                error("Error JHC09B in computeApproximateMDFordering.\n");
                            Int pikSpP = pijSpP;
                            double Cik = sys.C[pikSpP];
                            for (Int j = 1; j<nrk; j++) {
                                Int pkjSpP = ent2entStartP[rk]+j;
                                Int rj = ent2entP[pkjSpP];
                                if (j != i && sys.w[rj] < bign) {
                                    double Ckj = sys.C[pkjSpP];
                                    
                                    Int shared = 0;
                                    for (Int m = 1; m < nri; m++) {     // Determine if j neighbors i
                                        Int pimSpP = ent2entStartP[ri] + m;
                                        if (ent2entP[pimSpP] == rj) {
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
        for (Int rk=0; rk<nrowsP; rk++)
            sys.unordered2ordered[sys.ordered2unordered[rk]] = rk;
    
    delete[] pl_local;
}

void computeExactMDFordering(sysstruct &sys)
{
    /* TODO: Further optimization: If necessary, this function can be optimized by not recomputing w[i] when needs to be updated (we can just remove the contributions that no longer exist). */

    /* Note:
     * In this function, entity refers to faces in HDG and to face nodes in EDG and IEDG.
     * We assume all blocks are of the same size (i.e., nch (EDG) or nch*npf (UDH) are constant throughout the mesh).
     *      If not, an array indicating the size of each block and an array indicating the starting position each block in the matrix J would be required (easy fix).
     * The complexity of the MDF ordering is O(N^2) due to the indexofSmallestElement operation.
     * */
    
    Int inc = 1, info;
    Int nrowsP = sys.numEntitiesP;
    Int nblksP = sys.numBlocksP;
    Int *entSz = &sys.entSz[0];
    Int *ent2entStartP = &sys.ent2entStartP[0];
    Int *ent2entP = &sys.ent2entP[0];
    Int *rowStartP = &sys.rowStartP[0];
    Int *blkStartP = &sys.blkStartP[0];
    double one = 1.0, zero = 0.0;
    double bign = 1.0e100;
    char chn = 'N';

    vector<vector<Int> >* ipivs = &sys.ipivs;
    vector<vector<double> >* works = &sys.works;
    Int *pl_local = new Int[sys.noThreads];
    
    copyHg2MgWithMgFormat(sys);
    
//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for //private(i, j, m, ri, rj, nri, nrk, ik, kj, kk, pij, pim, pki, dnr, shared, info)
        for (Int rk=0; rk<nrowsP; rk++) {
            sys.w[rk] = 0.0;
            Int pkSpP = ent2entStart[rk];
            Int nrk = ent2entStartP[rk+1] - ent2entStartP[rk];
            Int pkP = rowStartP[rk];
            Int bsz_k = entSz[rk];
            Int bsz2 = bsz_k * bsz_k;
            Int lwork = bsz_k;
            
            // HkkInv <- inv(Hg_kk)
            DCOPY(&bsz2, &sys.Mg[pkP], &inc, &sys.HkkInvs[this_thread][0], &inc);
            DGETRF(&bsz_k, &bsz_k, &sys.HkkInvs[this_thread][0], &bsz_k, &ipivs[0][this_thread][0], &info);
            DGETRI(&bsz_k, &sys.HkkInvs[this_thread][0], &bsz_k, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
            for (Int i = 1; i<nrk; i++) {
                Int pkiSpP = pkSpP+i;
                Int ri = ent2entP[pkiSpP];
                Int nri = ent2entStartP[ri+1] - ent2entStartP[ri];
                Int bsz_i = entSz[ri];
                bsz2 = bsz_i * bsz_k;
                
                Int pijSpP;
                for (j=1; j<nri; j++) {
                    pijSpP = ent2entStartP[ri]+j;
                    if (ent2entP[pijSpP] == rk) break;
                }
                if (j == nri)
                    error("Error BS0DF9 in computeExactMDFordering.\n");
                Int pikSpP = pijSpP;
                Int pikP = blkStartP[pikSpP];
                
                /* HikHkkInv <- H_ik*HkkInv (=H_ik*inv(H_kk)) */
                DGEMM(&chn, &chn, &bsz_i, &bsz_k, &bsz_k, &one, &sys.Mg[pikP], &bsz_i, &sys.HkkInvs[this_thread][0],
                      &bsz_k, &zero, &sys.HikHkkInvs[this_thread][0], &bsz_i);
                for (Int j = 1; j<nrk; j++)
                    if (j != i) {
                        Int pkjSpP = pkSpP+j;
                        Int rj = ent2entP[pkjSpP];
                        Int bsz_j = entSz[rj];

                        Int shared = 0;
                        for (Int m=1; m<nri; m++) {     // Determine if j neighbors i
                            Int pimSpP = ent2entStartP[ri]+m;
                            if (ent2entP[pimSpP] == rj) {
                                shared = 1;
                                break;
                            }
                        }

                        if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                            Int pkjSpP = pkSpP + j;
                            Int pkjP = blkStartP[pkjSpP];
                            bsz2 = bsz_i * bsz_j;
                            
                            /* HikHkkInvHkj <- HikHkkInv*H_kj (=H_ik*inv(H_kk)*H_kj) */
                            DGEMM(&chn, &chn, &bsz_i, &bsz_j, &bsz_k, &one, &sys.HikHkkInvs[this_thread][0], &bsz_i, &sys.Mg[pkjP],
                                  &bsz_k, &zero, &sys.HikHkkInvHkjs[this_thread][0], &bsz_i);
                            double dnr = DNRM2(&bsz2, &sys.HikHkkInvHkjs[this_thread][0], &inc);
                            sys.w[rk] += dnr*dnr;   // TODO: Is this line correct or are we doing ^4?? (same with algorithm with constraints)
                        }
                    }
            }
            sys.w[rk] = sqrt(sys.w[rk]);
        }
        
        // Check NaN and Inf in sys.w:
        checkw(&sys.w[0], nrowsP, bign);
        
        for (Int rl=0; rl<nrowsP; rl++) {
// // //             Int chunkLen, startPt;
// // //             splitArray1D(nrowsP, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //             pl_local[this_thread] = startPt + indexofSmallestElement(&sys.w[startPt], chunkLen);
            Int pl = indexofSmallestElement(&sys.w[0], nrowsP);
// // //             #pragma omp barrier
// // //             #pragma omp single
// // //             {
// // //                 Int pl = pl_local[0];
// // //                 for (Int i_th = 1; i_th < sys.noThreads; i_th++)
// // //                     if (sys.w[pl_local[i_th]] < sys.w[pl])
// // //                         pl = pl_local[i_th];
                sys.ordered2unordered[rl] = pl;
                sys.w[pl] = bign;
                Int nrl = ent2entStartP[pl+1] - ent2entStartP[pl];
// // //             }
//             #pragma omp for //private(i, j, k, m, kk, ri, rj, rk, nri, nrk, pij, pki, plk, pim, ik, kj, dnr, shared, info)
            for (Int k = 1; k<nrl; k++) {
                Int plkSpP = ent2entStartP[pl]+k;
                Int rk = ent2entP[plkSpP];
                Int pkSpP = ent2entStartP[rk];
                Int pkP = rowStartP[rk];
                if (sys.w[rk] < bign) {
                    sys.w[rk] = 0.0;
                    Int nrk = ent2entStartP[rk+1] - ent2entStartP[rk];
                    Int bsz_k = entSz[rk];
                    Int bsz2 = bsz_k * bsz_k;
                    Int lwork = bsz_k;
                    
                    // HkkInv <- inv(Hg_kk)
                    DCOPY(&bsz2, &sys.Mg[pkP], &inc, &sys.HkkInvs[this_thread][0], &inc);
                    DGETRF(&bsz_k, &bsz_k, &sys.HkkInvs[this_thread][0], &bsz_k, &ipivs[0][this_thread][0], &info);
                    DGETRI(&bsz_k, &sys.HkkInvs[this_thread][0], &bsz_k, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
                    for (Int i = 1; i<nrk; i++) {
                        Int pkiSpP = pkSpP+i;
                        Int ri = ent2entP[pkiSpP];
                        if (sys.w[ri] < bign) {
                            Int bsz_i = entSz[ri];
                            Int nri = ent2entStartP[ri+1] - ent2entStartP[ri];
                            bsz2 = bsz_i * bsz_k;
                            
                            Int pijSpP;
                            for (Int j = 1; j< nri; j++) {
                                pijSpP = ent2entStartP[ri]+j;
                                if (ent2entP[pijSpP] == rk) break;
                            }
                            if (j == nri)
                                error("Error NJOI0L in computeExactMDFordering.\n");
                            Int pikSpP = pijSpP;
                            Int pikP = blkStartP[pikSpP];
                            
                            DGEMM(&chn, &chn, &bsz_i, &bsz_k, &bsz_k, &one, &sys.Mg[pikP], &bsz_i, &sys.HkkInvs[this_thread][0],
                                  &bsz_k, &zero, &sys.HikHkkInvs[this_thread][0], &bsz_i);
                            
                            for (Int j = 1; j<nrk; j++) {
                                Int pkjSpP = pkSpP+j;
                                Int rj = ent2entP[pkjSpP];
                                Int bsz_j = entSz[rj];
                                
                                if (j != i && sys.w[rj] < bign) {
                                    Int shared = 0;
                                    for (Int m = 1; m < nri; m++) {     // Determine if j neighbors i
                                        Int pimSpP = ent2entStartP[ri] + m;
                                        if (ent2entP[pimSpP] == rj) {
                                            shared = 1;
                                            break;
                                        }
                                    }
                                    if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                                        Int pkjSpP = pkSpP + j;
                                        Int pkjP = blkStartP[pkjSpP];
                                        bsz2 = bsz_i * bsz_j;
                                        
                                        DGEMM(&chn, &chn, &bsz_i, &bsz_j, &bsz_k, &one, &sys.HikHkkInvs[this_thread][0], &bsz_i, &sys.Mg[pkjP],
                                              &bsz_k, &zero, &sys.HikHkkInvHkjs[this_thread][0], &bsz_i);
                                        double dnr = DNRM2(&bsz2, &sys.HikHkkInvHkjs[this_thread][0], &inc);
                                        sys.w[rk] += dnr*dnr;   // TODO: Is this line correct or are we doing ^4?? (same with algorithm with constraints)
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
        for (Int rk=0; rk<nrowsP; rk++)
            sys.unordered2ordered[sys.ordered2unordered[rk]] = rk;
//     }
    
    delete[] pl_local;
}

void computeApproximateMDForderingWithConstraints(sysstruct &sys)
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
        error("\nApproximate MDF ordering with constraints only valid for parallel preconditioner.\n");
    }

    Int inc = 1, info;
    Int nrowsP = sys.numEntitiesP;
    Int nrowsPBJ = sys.BJ_nrowsP;
    Int nblksP = sys.numBlocksP;
    Int *entSz = &sys.entSz[0];
    Int *ent2entStartP = &sys.ent2entStartP[0];
    Int *ent2entP = &sys.ent2entP[0];
    Int *rowStartP = &sys.rowStartP[0];
    double one = 1.0, zero = 0.0;
    double bign = 1.0e100;
    char chn = 'N';

    vector<vector<Int> >* ipivs = &sys.ipivs;
    vector<vector<double> >* works = &sys.works;
    Int *pl_local = new Int[sys.noThreads];
    
    //copyHg2MgWithHgPlusFillFormat(sys); ??
    copyHg2MgWithMgFormat(sys); ??
    
//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for //private(info, j, ii, nri, pij)
        for (Int ri=0; ri<nrowsP; ri++) {
            Int piSpP = ent2entStartP[ri];
            Int pijP = rowStartP[ri];
            Int nri = ent2entStartP[ri+1] - ent2entStartP[ri];
            Int bsz_row = entSz[ri];
            Int bsz2 = bsz_row * bsz_row;
            Int lwork = bsz_row;
            
            /* HiiInv <- inv(Hg_ii) */
            DCOPY(&bsz2, &sys.Mg[pijP], &inc, &sys.HiiInvs[this_thread][0], &inc);
            DGETRF(&bsz_row, &bsz_row, &sys.HiiInvs[this_thread][0], &bsz_row, &ipivs[0][this_thread][0], &info);
            DGETRI(&bsz_row, &sys.HiiInvs[this_thread][0], &bsz_row, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
            
            for (Int j=0; j<nri; j++) {
                Int pijSpP = piSpP+j;
                Int rj = ent2entP[pijSpP];
                Int bsz_col = entSz[rj];
                Int bsz2 = bsz_row * bsz_col;
                
                /* HiiInvHij <- HiiInv*Hg_ij (=inv(Hg_ii)*Hg_ij)  */
                if (matchSparsity(ri,j)) {
                    DGEMM(&chn, &chn, &bsz_row, &bsz_col, &bsz_row, &one, &sys.HiiInvs[this_thread][0], &bsz_row, &sys.Mg[pijP],
                          &bsz_col, &zero, &sys.HiiInvHijs[this_thread][0], &bsz_row);
                    sys.C[pijSpP] = DNRM2(&bsz2, &sys.HiiInvHijs[this_thread][0], &inc);
                }
                else {
                    sys.C[pijSpP] = 0.0;
                }
                pijP += bsz2;
            }
        }
        
//         #pragma omp for //private(i, j, m, ri, rj, nri, nrk, ik, kj, kk, pij, pim, pki, Cik, Ckj, shared)
        for (Int rk=0; rk<nrowsP; rk++) {
            sys.w[rk] = 0.0;
            Int pkSpP = ent2entStartP[rk];
            Int pkiP = rowStartP[rk];
            Int nrk = ent2entStartP[rk+1] - ent2entStartP[rk];
            for (Int i = 1; i<nrk; i++) {
                Int pkiSpP = pkSpP + i;
                Int ri = ent2entP[pkiSpP];
                Int nri = ent2entStartP[ri+1] - ent2entStartP[ri];
                Int pijSpP;
                for (Int j=1; j<nri; j++) {
                    pijSpP = ent2entStartP[ri]+j;
                    if (ent2entP[pijSpP] == rk) break;
                }
                if (j == nri)
                    error("Error MXO9LJ in computeApproximateMDForderingWithConstraints.\n");
                Int pikSpP = pijSpP;
                double Cik = sys.C[pikSpP];
                for (Int j = 1; j<nrk; j++)
                    if (j != i) {
                        Int pkjSpP = pkSpP+j;
                        Int rj = ent2entP[pkjSpP];
                        double Ckj = sys.C[pkjSpP];

                        Int shared = 0;
                        for (Int m=1; m<nri; m++) {     // Determine if j neighbors i
                            Int pimSpP = ent2entStartP[ri]+m;
                            if (ent2entP[pimSpP] == rj) {
                                shared = 1;
                                break;
                            }
                        }
                        if (shared == 0)          // If j does not neighbor i, then increase the discarded fill in
                            sys.w[rk] += Cik*Cik*Ckj*Ckj;
                    }
            }
            sys.w[rk] = sqrt(sys.w[rk]);
        }
//     }
        
        // Check NaN and Inf in sys.w:
        checkw(&sys.w[0], nrowsP, bign);
        
        for (Int rl=0; rl<nrowsP; rl++) {
// // //             Int chunkLen, startPt;
// // //             if (rl < nrowsPBJ) {
// // //                 splitArray1D(nrowsPBJ, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //                 pl_local[this_thread] = startPt + indexofSmallestElement(&sys.w[startPt], chunkLen);
// // //             }
// // //             else {
// // //                 splitArray1D(nrowsP-nrowsPBJ, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //                 pl_local[this_thread] = nrowsPBJ+startPt + indexofSmallestElement(&sys.w[nrowsPBJ+startPt], chunkLen);
// // //             }
            Int pl;
            if (rl < nrowsPBJ)
                pl = indexofSmallestElement(&sys.w[0], nrowsPBJ);
            else
                pl = nrowsPBJ + indexofSmallestElement(&sys.w[nrowsPBJ], nrowsP-nrowsPBJ);
// // //             #pragma omp barrier
// // //             #pragma omp single
// // //             {
// // //                 pl = pl_local[0];
// // //                 for (Int i_th = 1; i_th < sys.noThreads; i_th++)
// // //                     if (sys.w[pl_local[i_th]] < sys.w[pl])
// // //                         pl = pl_local[i_th];
                sys.ordered2unordered[rl] = pl;
                sys.w[pl] = bign;
                Int nrl = ent2entStartP[pl+1] - ent2entStartP[pl];
// // //             }
// // //             #pragma omp for //private(i, j, k, m, ri, rj, rk, nri, nrk, pij, pki, plk, pim, ik, kj, Cik, Ckj, shared)
            for (Int k = 1; k<nrl; k++) {
                Int plkSpP = ent2entStartP[pl]+k;
                Int rk = ent2entP[plkSpP];
                if (sys.w[rk] < bign) {
                    sys.w[rk] = 0.0;
                    Int nrk = ent2entStartP[rk+1] - ent2entStartP[rk];
                    for (Int i = 1; i<nrk; i++) {
                        Int pkiSpP = ent2entStartP[rk]+i;
                        Int ri = ent2entP[pkiSpP];
                        if (sys.w[ri] < bign) {
                            Int nri = ent2entStartP[ri+1] - ent2entStartP[ri];
                            Int pijSpP;
                            for (Int j = 1; j<nri; j++) {
                                pijSpP = ent2entStartP[ri]+j;
                                if (ent2entP[pijSpP] == rk) break;
                            }
                            if (j == nri)
                                error("Error E89J09 in computeApproximateMDForderingWithConstraints.\n");
                            Int pikSpP = pijSpP;
                            double Cik = sys.C[pikSpP];
                            for (Int j = 1; j<nrk; j++) {
                                Int pkjSpP = ent2entStartP[rk]+j;
                                Int rj = ent2entP[pkjSpP];
                                if (j != i && sys.w[rj] < bign) {
                                    double Ckj = sys.C[pkjSpP];
                                    
                                    Int shared = 0;
                                    for (Int m = 1; m < nri; m++) {     // Determine if j neighbors i
                                        Int pimSpP = ent2entStartP[ri] + m;
                                        if (ent2entP[pimSpP] == rj) {
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
        for (Int rk=0; rk<nrowsP; rk++)
            sys.unordered2ordered[sys.ordered2unordered[rk]] = rk;
    
    delete[] pl_local;
}

void computeExactMDForderingWithConstraints(sysstruct &sys)
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
        error("\nExact MDF ordering with constraints only valid for parallel preconditioner.\n");
    }
    
    Int inc = 1, info;
    Int nrowsP = sys.numEntitiesP;
    Int nrowsPBJ = sys.BJ_nrowsP;
    Int nblksP = sys.numBlocksP;
    Int *entSz = &sys.entSz[0];
    Int *ent2entStartP = &sys.ent2entStartP[0];
    Int *ent2entP = &sys.ent2entP[0];
    Int *rowStartP = &sys.rowStartP[0];
    Int *blkStartP = &sys.blkStartP[0];
    double one = 1.0, zero = 0.0;
    double bign = 1.0e100;
    char chn = 'N';

    vector<vector<Int> >* ipivs = &sys.ipivs;
    vector<vector<double> >* works = &sys.works;
    Int *pl_local = new Int[sys.noThreads];
    
    copyHg2MgWithMgFormat(sys);
    
//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for //private(i, j, m, ri, rj, nri, nrk, ik, kj, kk, pij, pim, pki, dnr, shared, info)
        for (Int rk=0; rk<nrowsP; rk++) {
            sys.w[rk] = 0.0;
            Int pkSpP = ent2entStart[rk];
            Int nrk = ent2entStartP[rk+1] - ent2entStartP[rk];
            Int pkP = rowStartP[rk];
            Int bsz_k = entSz[rk];
            Int bsz2 = bsz_k * bsz_k;
            Int lwork = bsz_k;
            
            // HkkInv <- inv(Hg_kk)
            DCOPY(&bsz2, &sys.Mg[pkP], &inc, &sys.HkkInvs[this_thread][0], &inc);
            DGETRF(&bsz_k, &bsz_k, &sys.HkkInvs[this_thread][0], &bsz_k, &ipivs[0][this_thread][0], &info);
            DGETRI(&bsz_k, &sys.HkkInvs[this_thread][0], &bsz_k, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
            for (Int i = 1; i<nrk; i++) {
                Int pkiSpP = pkSpP+i;
                Int ri = ent2entP[pkiSpP];
                Int nri = ent2entStartP[ri+1] - ent2entStartP[ri];
                Int bsz_i = entSz[ri];
                bsz2 = bsz_i * bsz_k;
                
                Int pijSpP;
                for (j=1; j<nri; j++) {
                    pijSpP = ent2entStartP[ri]+j;
                    if (ent2entP[pijSpP] == rk) break;
                }
                if (j == nri)
                    error("Error M3207H in computeExactMDForderingWithConstraints.\n");
                Int pikSpP = pijSpP;
                Int pikP = blkStartP[pikSpP];
                
                /* HikHkkInv <- H_ik*HkkInv (=H_ik*inv(H_kk)) */
                DGEMM(&chn, &chn, &bsz_i, &bsz_k, &bsz_k, &one, &sys.Mg[pikP], &bsz_i, &sys.HkkInvs[this_thread][0],
                      &bsz_k, &zero, &sys.HikHkkInvs[this_thread][0], &bsz_i);
                for (Int j = 1; j<nrk; j++)
                    if (j != i) {
                        Int pkjSpP = pkSpP+j;
                        Int rj = ent2entP[pkjSpP];
                        Int bsz_j = entSz[rj];

                        Int shared = 0;
                        for (Int m=1; m<nri; m++) {     // Determine if j neighbors i
                            Int pimSpP = ent2entStartP[ri]+m;
                            if (ent2entP[pimSpP] == rj) {
                                shared = 1;
                                break;
                            }
                        }

                        if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                            Int pkjSpP = pkSpP + j;
                            Int pkjP = blkStartP[pkjSpP];
                            bsz2 = bsz_i * bsz_j;
                            
                            /* HikHkkInvHkj <- HikHkkInv*H_kj (=H_ik*inv(H_kk)*H_kj) */
                            DGEMM(&chn, &chn, &bsz_i, &bsz_j, &bsz_k, &one, &sys.HikHkkInvs[this_thread][0], &bsz_i, &sys.Mg[pkjP],
                                  &bsz_k, &zero, &sys.HikHkkInvHkjs[this_thread][0], &bsz_i);
                            double dnr = DNRM2(&bsz2, &sys.HikHkkInvHkjs[this_thread][0], &inc);
                            sys.w[rk] += dnr*dnr;   // TODO: Is this line correct or are we doing ^4?? (same with algorithm with constraints)
                        }
                    }
            }
            sys.w[rk] = sqrt(sys.w[rk]);
        }
        
        // Check NaN and Inf in sys.w:
        checkw(&sys.w[0], nrowsP, bign);
        
        for (Int rl=0; rl<nrowsP; rl++) {
// // //             Int chunkLen, startPt;
// // //             if (rl < nrowsPBJ) {
// // //                 splitArray1D(nrowsPBJ, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //                 pl_local[this_thread] = startPt + indexofSmallestElement(&sys.w[startPt], chunkLen);
// // //             }
// // //             else {
// // //                 splitArray1D(nrowsP-nrowsPBJ, (Int) this_thread, sys.noThreads, &chunkLen, &startPt);
// // //                 pl_local[this_thread] = nrowsPBJ+startPt + indexofSmallestElement(&sys.w[nrowsPBJ+startPt], chunkLen);
// // //             }
            Int pl;
            if (rl < nrowsPBJ)
                pl = indexofSmallestElement(&sys.w[0], nrowsPBJ);
            else
                pl = nrowsPBJ + indexofSmallestElement(&sys.w[nrowsPBJ], nrowsP-nrowsPBJ);
// // //             #pragma omp barrier
// // //             #pragma omp single
// // //             {
// // //                 pl = pl_local[0];
// // //                 for (Int i_th = 1; i_th < sys.noThreads; i_th++)
// // //                     if (sys.w[pl_local[i_th]] < sys.w[pl])
// // //                         pl = pl_local[i_th];
                sys.ordered2unordered[rl] = pl;
                sys.w[pl] = bign;
                Int nrl = ent2entStartP[pl+1] - ent2entStartP[pl];
// // //             }
//             #pragma omp for //private(i, j, k, m, kk, ri, rj, rk, nri, nrk, pij, pki, plk, pim, ik, kj, dnr, shared, info)
            for (Int k = 1; k<nrl; k++) {
                Int plkSpP = ent2entStartP[pl]+k;
                Int rk = ent2entP[plkSpP];
                Int pkSpP = ent2entStartP[rk];
                Int pkP = rowStartP[rk];
                if (sys.w[rk] < bign) {
                    sys.w[rk] = 0.0;
                    Int nrk = ent2entStartP[rk+1] - ent2entStartP[rk];
                    Int bsz_k = entSz[rk];
                    Int bsz2 = bsz_k * bsz_k;
                    Int lwork = bsz_k;
                    
                    // HkkInv <- inv(Hg_kk)
                    DCOPY(&bsz2, &sys.Mg[pkP], &inc, &sys.HkkInvs[this_thread][0], &inc);
                    DGETRF(&bsz_k, &bsz_k, &sys.HkkInvs[this_thread][0], &bsz_k, &ipivs[0][this_thread][0], &info);
                    DGETRI(&bsz_k, &sys.HkkInvs[this_thread][0], &bsz_k, &ipivs[0][this_thread][0], &works[0][this_thread][0], &lwork, &info);
                    for (Int i = 1; i<nrk; i++) {
                        Int pkiSpP = pkSpP+i;
                        Int ri = ent2entP[pkiSpP];
                        if (sys.w[ri] < bign) {
                            Int bsz_i = entSz[ri];
                            Int nri = ent2entStartP[ri+1] - ent2entStartP[ri];
                            bsz2 = bsz_i * bsz_k;
                            
                            Int pijSpP;
                            for (Int j = 1; j< nri; j++) {
                                pijSpP = ent2entStartP[ri]+j;
                                if (ent2entP[pijSpP] == rk) break;
                            }
                            if (j == nri)
                                error("Error PW8G7K in computeExactMDForderingWithConstraints.\n");
                            Int pikSpP = pijSpP;
                            Int pikP = blkStartP[pikSpP];
                            
                            DGEMM(&chn, &chn, &bsz_i, &bsz_k, &bsz_k, &one, &sys.Mg[pikP], &bsz_i, &sys.HkkInvs[this_thread][0],
                                  &bsz_k, &zero, &sys.HikHkkInvs[this_thread][0], &bsz_i);
                            
                            for (Int j = 1; j<nrk; j++) {
                                Int pkjSpP = pkSpP+j;
                                Int rj = ent2entP[pkjSpP];
                                Int bsz_j = entSz[rj];
                                
                                if (j != i && sys.w[rj] < bign) {
                                    Int shared = 0;
                                    for (Int m = 1; m < nri; m++) {     // Determine if j neighbors i
                                        Int pimSpP = ent2entStartP[ri] + m;
                                        if (ent2entP[pimSpP] == rj) {
                                            shared = 1;
                                            break;
                                        }
                                    }
                                    if (shared == 0) {          // If j does not neighbor i, then increase the discarded fill in
                                        Int pkjSpP = pkSpP + j;
                                        Int pkjP = blkStartP[pkjSpP];
                                        bsz2 = bsz_i * bsz_j;
                                        
                                        DGEMM(&chn, &chn, &bsz_i, &bsz_j, &bsz_k, &one, &sys.HikHkkInvs[this_thread][0], &bsz_i, &sys.Mg[pkjP],
                                              &bsz_k, &zero, &sys.HikHkkInvHkjs[this_thread][0], &bsz_i);
                                        double dnr = DNRM2(&bsz2, &sys.HikHkkInvHkjs[this_thread][0], &inc);
                                        sys.w[rk] += dnr*dnr;   // TODO: Is this line correct or are we doing ^4?? (same with algorithm with constraints)
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
        for (Int rk=0; rk<nrowsP; rk++)
            sys.unordered2ordered[sys.ordered2unordered[rk]] = rk;
//     }
    
    delete[] pl_local;
}

void postprocessOrderingForBILUk(sysstruct &sys)
{
    // TODO: This function is hard to thread and the benefit may be small.
    
    Int nrowsP = sys.numRowsP;
    
    Int *entSz = &sys.entSz[0];
    Int *ent2entStartP = &sys.ent2entStartP[0];
    Int *ent2entP = &sys.ent2entP[0];
    Int *oent2entStart = &sys.oent2entStart[0];
    Int *oent2ent = &sys.oent2ent[0];
    Int *ordered2unordered = &sys.ordered2unordered[0];
    Int *unordered2ordered = &sys.unordered2ordered[0];
    Int *Loent2entStart = &sys.Loent2entStart[0];
    Int *Uoent2entStart = &sys.Uoent2entStart[0];
    Int *LUoent2ent = &sys.LUoent2ent[0];
    Int *UorowStart = &sys.UorowStart[0];
    Int *LorowStart = &sys.LorowStart[0];
    Int *LUoblkStart = &sys.LUoblkStart[0];
    
    // Compute sys.oent2entStart and sys.oent2ent:
    oent2entStart[0] = 0;
    for (Int rj=0; rj<nrowsP; rj++) {
        Int rjj = ordered2unordered[rj];
        Int pj = oent2entStart[rj];
        Int pjj = ent2entStartP[rjj];
        Int nrj = ent2entStartP[rjj+1]-ent2entStartP[rjj];
        
        oent2entStart[rj+1] = oent2entStart[rj] + nrj;
        for (Int i=0; i<nrj; i++)
            oent2ent[pj+i] = ent2entP[pjj+i];
    }
    
    // Compute sys.LUoent2ent, sys.LUoblkStart (parts associated to L), sys.Loent2entStart and sys.LorowStart:
    Loent2entStart[0] = 0;
    LorowStart[0] = 0;
    LUoblkStart[0] = 0;
    for (Int rj=0; rj<nrowsP; rj++) {
        Int rjj = ordered2unordered[rj];
        Int pj_l = Loent2entStart[rj];
        Int pjj = ent2entStartP[rjj];
        Int nrj = ent2entStartP[rjj+1]-ent2entStartP[rjj];
        Int bsz_row = entSz[rjj];
        
        Int i_l = 0, bszs = 0;
        for (Int i=1; i<nrj; i++) {
            Int pjji = pjj+i;
            Int rii = ent2entP[pjji];
            Int ri = unordered2ordered[rii];

            if (ri < rj) {
                Int bsz_col = entSz[rii];
                Int bsz2 = bsz_row * bsz_col;
                LUoent2ent[pj_l+i_l] = rii;
                LUoblkStart[pj_l+i_l+1] = LUoblkStart[pj_l+i_l] + bsz2;
                i_l ++;
                bszs += bsz2;
                
            }
        }
        Loent2entStart[rj+1] = Loent2entStart[rj] + i_l;
        LorowStart[rj+1] = LorowStart[rj] + bszs;
    }
    
    // Compute sys.Uoent2entStart and sys.UorowStart
    Uoent2entStart[nrowsP] = Loent2entStart[nrowsP];
    UorowStart[nrowsP] = LorowStart[nrowsP];
    for (Int j=0; j<nrowsP; j++) {
        Int rj = nrowsP-j-1;
        Int rjj = ordered2unordered[rj];
        Int pjj = ent2entStartP[rjj];
        Int nrj = ent2entStartP[rjj+1]-ent2entStartP[rjj];
        Int bsz_row = entSz[rjj];
        
        Int i_u;
        if (nrj > 0)
            i_u = 1;        // Start at 1 to consider the diagonal block
        else {              // This is a ghost (missing) entity in the mesh
            printf("\nWARNING: Ghost entities have been detected in the mesh.\n\n");
            i_u = 0;
        }
        
        Int bszs = 0;
        for (Int i=1; i<nrj; i++) {
            Int pjji = pjj+i;
            Int rii = ent2entP[pjji];
            Int ri = unordered2ordered[rii];
            if (ri > rj) {
                Int bsz_col = entSz[rii];
                i_u ++;
                bszs += bsz_row * bsz_col;
            }
        }
        Uoent2entStart[rj] = Uoent2entStart[rj+1] + i_u;
        UorowStart[rj] = UorowStart[rj+1] + bszs;
    }
    
    // Compute sys.LUoent2ent and sys.LUoblkStart (parts associated to U)
    for (Int j=0; j<nrowsP; j++) {
        Int rj = nrowsP-j-1;
        Int rjj = ordered2unordered[rj];
        Int pj_u = Uoent2entStart[rj+1];
        Int pjj = ent2entStartP[rjj];
        Int nrj = ent2entStartP[rjj+1] - ent2entStartP[rjj];
        Int bsz_row = entSz[rjj];

        Int i_u = 0;
        for (Int i=1; i<nrj; i++) {
            Int pjji = pjj+i;
            Int rii = ent2entP[pjji];
            Int ri = unordered2ordered[rii];
            if (ri > rj) {
                Int bsz_col = entSz[rii];
                Int bsz2 = bsz_row * bsz_col;
                LUoent2ent[pj_u+i_u] = rii;
                LUoblkStart[pj_u+i_u+1] = LUoblkStart[pj_u+i_u] + bsz2;
                i_u ++;
            }
        }
        LUoent2ent[pj_u+i_u] = rjj;
        LUoblkStart[pj_u+i_u+1] = LUoblkStart[pj_u+i_u] + bsz_row*bsz_row;
    }
}

void computeOrdering(sysstruct &sys)
{
    // TODO: Change the code so that for subdomainPreconditioner != 0 no reorder is even computed
    if (sys.subdomainPreconditioner != 0 && sys.reorderMethod != 0)
        sys.reorderMethod = 0;
    
    if (sys.reorderMethod == 0)             // No reorder
        computeNoOrdering(sys);
    else if (sys.reorderMethod == 1)        // Approximate MDF
        computeApproximateMDFordering(sys);
    else if (sys.reorderMethod == 2)        // Exact MDF
        computeExactMDFordering(sys);
    else if (sys.reorderMethod == 3)        // Approximate MDF with restricted order for entities not in subdomain. For RAS(l), l>0 parallel preconditioner only.
        computeApproximateMDForderingWithConstraints(sys);
    else if (sys.reorderMethod == 4)        // Exact MDF with restricted order for entities not in subdomain. For RAS(l), l>0 parallel preconditioner only.
        computeExactMDForderingWithConstraints(sys);
    else
        printf("Ordering algorithm not implemented yet. sys.reorderMethod = %d.\n",sys.reorderMethod);
    
    if (sys.subdomainPreconditioner == 0)
        postprocessOrderingForBILUk(sys);
}

#endif
