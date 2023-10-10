#ifndef __PSLAK
#define __PSLAK

// Written by: C. Nguyen & P. Fernandez

#ifdef  HAVE_MPI

// static void matvecfreeMPI(double* v, sysstruct &sys, double* x, Int sendInputVector)
// {
//     // TODO: Thread this function
//     
//     // Note 1: Same implementation both for single and double precision matrix-vector product
//     // Note 2: The current implementation is a bit hacky. udg and uh are used as 
//     //     temporary vectors to store sol.UDG and sol.UH since computeSchurResidualMPI
//     //     uses sol.UDG and sol.UH to compute the new residual.
//     // Note 3: sys.x may have been sent somewhere in the code and may need to be received here.
//
//     Int ne = ndims[5];
//     Int npv = ndims[9];
//     Int ncu = ndims[20];
//     Int n1 = npv*ncu;
//     
//     double *udg = &sol.UDG4FD[0];
//     double *uh = &sol.UH4FD[0];
//     double *duh = &x[0];
//     double *Rg1 = &sys.Rg4FD[0];
//     double *r1 = &sys.r4FD[0];
//     double *Rg0 = &sys.Rg_0[0];
//     double epsilon = 1.0e-6;
//     
//     for (Int i = 0; i < sol.UDG.size(); i++)
//         udg[i] = sol.UDG[i];
//     for (Int i = 0; i < sol.UH.size(); i++)
//         uh[i] = sol.UH[i];
//     
//     if (sendInputVector == 0)
//         updateSol(sys, elems, mesh, master, sol, app, temps, ndims, &sol.UDG[0], &sol.UH[0], duh, epsilon);
//     else if (sendInputVector == 1)
//         updateSolMPI(sys, elems, mesh, master, sol, app, temps, ndims, &sol.UDG[0], &sol.UH[0], duh, epsilon);
// 
//     computeSchurResidualMPI(sys, elems, mesh, master, sol, app, temps, ndims, &Rg1[0], &r1[0]);
//     
//     for (Int i = 0; i < sol.UDG.size(); i++)
//         sol.UDG[i] = udg[i];
//     for (Int i = 0; i < sol.UH.size(); i++)
//         sol.UH[i] = uh[i];
//     
//     Int n1 = sys.blkSize*sys.BJ_nrows;
//     for (Int i = 0; i < n1; i++)
//         v[i] = (Rg1[i] - Rg0[i]) / epsilon;
// }

// static void matvecMPI_BJ_sp(double* v, sysstruct &sys, double* x, Int * requestCounter, double* matvectimes)
// {
//     // Description: Single-precision, parallel matrix-vector product 
//     //              for block Jacobi -RAS(0)- preconditioner
// 
//     Int inc = 1, i, j, k, l, n, ri, rj, rm, nb, nri, pij, pii;
//     Int bsz = sys.blkSize;
//     Int bsz2 = bsz*bsz;
//     Int * BJ_rowpts = &sys.BJ_rowpts[0];
//     Int * BK_rowpts = &sys.BK_rowpts[0];
//     Int * BK_rowind = &sys.BK_rowind[0];
//     Int * BK_colind = &sys.BK_colind[0];
//     Int * BJ_colind = &sys.BJ_colind[0];
//     float * Hg_sp = &sys.Hg_sp[0];
//     float * x_sp = &sys.x_sp[0];
//     float * v_sp = &sys.v_sp[0];
//     float one = 1.0, zero = 0.0, alpha;
//     double matvectime1, matvectime2, MPIwaitTime, tTotal;
//     char chn = 'N';
//     clock_t t, t0 = clock();
//     
//     /* Convert x to single precision */
//     for (i = 0; i < sys.x_sp.size(); i++)
//         x_sp[i] = (float) x[i];
//     
//     /* perform v <- Hg*x while sending and receiving data */
//     t = clock();
//     //pij = -1;
//     if (sys.matvecImplementation == 0) {
// //         #pragma omp parallel for num_threads(sys.noThreads) private(alpha, nri, rj, pij, i)
//         for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
//             alpha = 0.0;
//             nri = BJ_rowpts[ri+1] - BJ_rowpts[ri];
// 
//             for (i = 0; i<nri; i++) {
//                 //pij += 1;
//                 pij = BJ_rowpts[ri]+i;
//                 rj = BJ_colind[pij];
// 
//                 SGEMV(&chn, &bsz, &bsz, &one, &Hg_sp[bsz2*pij], &bsz, &x_sp[bsz*rj],
//                       &inc, &alpha, &v_sp[bsz*ri], &inc);
// 
//                 alpha = 1.0;
//             }
//         }
//     }
//     else if (sys.matvecImplementation == 1) {
// //         #pragma omp parallel num_threads(sys.noThreads) private(nri, pii, pij, rj, nb, i, j)
//         {
// //             int this_thread = omp_get_thread_num();
//             int this_thread = 0;
// //             #pragma omp for
//             for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
//                 nri = BJ_rowpts[ri+1] - BJ_rowpts[ri];
//                 pii = BJ_rowpts[ri];
// 
//                 // Create temporary dense vector:
//                 for (i = 0; i<nri; i++) {
//                     //pij += 1;
//                     pij = BJ_rowpts[ri]+i;
//                     rj = BJ_colind[pij];
//                    for (j = 0; j < bsz; j++)
//                        sys.xDenseRow_sp[this_thread][bsz*i+j] = x_sp[bsz*rj+j];
//     //                 SCOPY(&bsz, &x_sp[bsz*rj], &inc, &sys.xDenseRow_sp[this_thread][bsz*i], &inc);
//                 }
// 
//                 // Perform matrix-vector product:
//                 nb = bsz*nri;
// 
//                 SGEMV(&chn, &bsz, &nb, &one, &Hg_sp[bsz2*pii], &bsz, &sys.xDenseRow_sp[this_thread][0],
//                       &inc, &zero, &v_sp[bsz*ri], &inc);
//             }
//         }
//     }
//     matvectimes[0] += clock() - t;
// 
//     /* Wait until all send and receive operations are completely done */
//     t = clock();
//     MPI_Waitall(*requestCounter, sys.requests, sys.statuses);
//     matvectimes[1] += clock() - t;
// 
//     /* Copy buffrecv to x_sp */
//     for (i=0; i<sys.nentrecv; i++) {
//         ri = sys.entrecv[i];
//         for (j=0; j<bsz; j++)
//             x_sp[bsz*ri+j] = (float) sys.buffrecv[bsz*i+j];
//     }
//     
//     /* v <- v + Kg*x */
//     t = clock();
//     //pij = -1;
//     if (sys.matvecImplementation == 0) {
// //         #pragma omp parallel for num_threads(sys.noThreads) private(nri, rj, rm, pij, i)
//         for (Int ri = 0; ri < sys.BK_nrows; ri++) {
//             rm = BK_rowind[ri];
//             nri = BK_rowpts[ri+1] - BK_rowpts[ri];
//             for (i = 0; i<nri; i++) {
//                 //pij += 1;
//                 pij = BK_rowpts[ri]+i;
//                 rj = BK_colind[pij];
//                 SGEMV(&chn, &bsz, &bsz, &one, &sys.Kg_sp[bsz2*pij], &bsz, &x_sp[bsz*rj],
//                       &inc, &one, &v_sp[bsz*rm], &inc);
//             }
//         }
//     }
//     else if (sys.matvecImplementation == 1) {
// //         #pragma omp parallel num_threads(sys.noThreads) private(nri, pii, pij, rj, rm, nb, i, j)
//         {
// //             int this_thread = omp_get_thread_num();
//             int this_thread = 0;
// //             #pragma omp for
//             for (Int ri = 0; ri < sys.BK_nrows; ri++) {
//                 rm = BK_rowind[ri];
//                 nri = BK_rowpts[ri+1] - BK_rowpts[ri];
//                 pii = BK_rowpts[ri];
// 
//                 // Create temporary dense vector:
//                 for (i = 0; i<nri; i++) {
//                     //pij += 1;
//                     pij = BK_rowpts[ri]+i;
//                     rj = BK_colind[pij];
//                    for (j = 0; j < bsz; j++)
//                        sys.xDenseRow_sp[this_thread][bsz*i+j] = x_sp[bsz*rj+j];
//     //                 SCOPY(&bsz, &x_sp[bsz*rj], &inc, &sys.xDenseRow_sp[this_thread][bsz*i], &inc);
//                 }
// 
//                 // Perform matrix-vector product:
//                 nb = bsz*nri;
//                 SGEMV(&chn, &bsz, &nb, &one, &sys.Kg_sp[bsz2*pii], &bsz, &sys.xDenseRow_sp[this_thread][0],
//                       &inc, &one, &v_sp[bsz*rm], &inc);
//             }
//         }
//     }
//     matvectimes[2] += clock() - t;
//     
//     /* Convert v to double precision */
//     for (i = 0; i < sys.v_sp.size(); i++)
//         v[i] = (double) v_sp[i];
//     
//     tTotal = clock() - t0;
// }
// 
// static void matvecMPI_BJ_dp(double* v, sysstruct &sys, double* x, Int * requestCounter, double* matvectimes)
// {
//     // Description: Double-precision, parallel matrix-vector product 
//     //              for block Jacobi -RAS(0)- preconditioner
// 
//     Int i, j, k, l, n, ri, rj, rm, inc = 1;
//     Int bsz = sys.blkSize;
//     Int bsz2 = bsz*bsz;
//     Int nb;
//     Int nri, pij, pii;
//     Int * BJ_rowpts = &sys.BJ_rowpts[0];
//     Int * BK_rowpts = &sys.BK_rowpts[0];
//     Int * BK_rowind = &sys.BK_rowind[0];
//     Int * BK_colind = &sys.BK_colind[0];
//     Int * BJ_colind = &sys.BJ_colind[0];
//     double * Hg = &sys.Hg[0];
//     double one = 1.0, zero = 0.0, alpha;
//     double matvectime1, matvectime2, MPIwaitTime, tTotal;
//     char chn = 'N';
//     clock_t t, t0 = clock();
//     
//     /* perform v <- Hg*x while sending and receiving data */
//     t = clock();
//     //pij = -1;
//     if (sys.matvecImplementation == 0) {
// //         #pragma omp parallel for num_threads(sys.noThreads) private(alpha, nri, rj, pij, i)
//         for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
//             alpha = 0.0;
//             nri = BJ_rowpts[ri+1] - BJ_rowpts[ri];
// 
//             for (i = 0; i<nri; i++) {
//                 //pij += 1;
//                 pij = BJ_rowpts[ri]+i;
//                 rj = BJ_colind[pij];
// 
//                 DGEMV(&chn, &bsz, &bsz, &one, &Hg[bsz2*pij], &bsz, &x[bsz*rj],
//                       &inc, &alpha, &v[bsz*ri], &inc);
// 
//                 alpha = 1.0;
//             }
//         }
//     }
//     else if (sys.matvecImplementation == 1) {
// //         #pragma omp parallel num_threads(sys.noThreads) private(nri, pii, pij, rj, nb, i, j)
//         {
// //             int this_thread = omp_get_thread_num();
//             int this_thread = 0;
// //             #pragma omp for
//             for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
//                 nri = BJ_rowpts[ri+1] - BJ_rowpts[ri];
//                 pii = BJ_rowpts[ri];
// 
//                 // Create temporary dense vector:
//                 for (i = 0; i<nri; i++) {
//                     //pij += 1;
//                     pij = BJ_rowpts[ri]+i;
//                     rj = BJ_colind[pij];
//                    for (j = 0; j < bsz; j++)
//                        sys.xDenseRow[this_thread][bsz*i+j] = x[bsz*rj+j];
//                     //DCOPY(&bsz, &x[bsz*rj], &inc, &sys.xDenseRow[this_thread][bsz*i], &inc);
//                 }
//                 
//                 // Perform matrix-vector product
//                 nb = bsz*nri;
// 
//                 DGEMV(&chn, &bsz, &nb, &one, &Hg[bsz2*pii], &bsz, &sys.xDenseRow[this_thread][0],
//                       &inc, &zero, &v[bsz*ri], &inc);
//             }
//         }
//     }
//     matvectimes[0] += clock() - t;
// 
//     /* wait until all send and receive operations are completely done */
//     t = clock();
//     MPI_Waitall(*requestCounter, sys.requests, sys.statuses);
//     matvectimes[1] += clock() - t;
// 
//     /* copy buffrecv to x */
//     for (i=0; i<sys.nentrecv; i++) {
//         ri = sys.entrecv[i];
//         DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &x[bsz*ri], &inc);
//     }
//     
//     /* v <- v + Kg*x */
//     t = clock();
//     //pij = -1;
//     if (sys.matvecImplementation == 0) {
// //         #pragma omp parallel for num_threads(sys.noThreads) private(nri, rj, rm, pij, i)
//         for (Int ri = 0; ri < sys.BK_nrows; ri++) {
//             rm = BK_rowind[ri];
//             nri = BK_rowpts[ri+1] - BK_rowpts[ri];
//             for (i = 0; i<nri; i++) {
//                 //pij += 1;
//                 pij = BK_rowpts[ri]+i;
//                 rj = BK_colind[pij];
//                 DGEMV(&chn, &bsz, &bsz, &one, &sys.Kg[bsz2*pij], &bsz, &x[bsz*rj],
//                       &inc, &one, &v[bsz*rm], &inc);
//             }
//         }
//     }
//     else if (sys.matvecImplementation == 1) {
// //         #pragma omp parallel num_threads(sys.noThreads) private(nri, pii, pij, rj, rm, nb, i, j)
//         {
// //             int this_thread = omp_get_thread_num();
//             int this_thread = 0;
// //             #pragma omp for
//             for (Int ri = 0; ri < sys.BK_nrows; ri++) {
//                 rm = BK_rowind[ri];
//                 nri = BK_rowpts[ri+1] - BK_rowpts[ri];
//                 pii = BK_rowpts[ri];
// 
//                 // Create temporary dense vector:
//                 for (i = 0; i<nri; i++) {
//                     //pij += 1;
//                     pij = BK_rowpts[ri]+i;
//                     rj = BK_colind[pij];
//                    for (j = 0; j < bsz; j++)
//                        sys.xDenseRow[this_thread][bsz*i+j] = x[bsz*rj+j];
//                     //DCOPY(&bsz, &x[bsz*rj], &inc, &sys.xDenseRow[this_thread][bsz*i], &inc);
//                 }
// 
//                 // Perform matrix-vector product
//                 nb = bsz*nri;
//                 DGEMV(&chn, &bsz, &nb, &one, &sys.Kg[bsz2*pii], &bsz, &sys.xDenseRow[this_thread][0],
//                       &inc, &one, &v[bsz*rm], &inc);
//             }
//         }
//     }
//     matvectimes[2] += clock() - t;
//     
//     tTotal = clock() - t0;
// }

static void matvecMPI_RAS_sp(double* v, sysstruct &sys, double* x, Int * requestCounter, double* matvectimes)
{
    // Description: Single-precision, parallel matrix-vector product 
    //              for RAS(1) preconditioner
    
    Int i, j, n, nri, nriInt, nriExt, pij, pii, pim, ri, rj, rij, inc = 1;
    Int nb;
    Int bsz = sys.blkSize;
    Int bsz2    = bsz*bsz;
    Int *ent2entStart = &sys.ent2entStart[0];
    Int *ent2ent = &sys.ent2ent[0];
    float *Hg_sp = &sys.Hg_sp[0];
    float *x_sp = &sys.x_sp[0];
    float *v_sp = &sys.v_sp[0];
    double matvectime1, matvectime2, MPIwaitTime, tTotal;
    float one = 1.0, zero = 0.0, alpha;
    char chn = 'N';
    clock_t t, t0 = clock();

    /* Convert x to single precision */
    for (i = 0; i < sys.x_sp.size(); i++)
        x_sp[i] = (float) x[i];

    /* v <- Hg*x for entities in the processor while sending and receiving data */
    t = clock();

    switch (sys.matvecImplementation) {
    case 0:
//         #pragma omp parallel for num_threads(sys.noThreads) private(alpha, nri, rj, pij, i)
        for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
            alpha = 0.0;
            nri = ent2entStart[ri+1] - ent2entStart[ri];
            for (i = 0; i<nri; i++) {
                pij = ent2entStart[ri]+i;
                rj = ent2ent[pij];
                if (rj < sys.BJ_nrows) {
                    SGEMV(&chn, &bsz, &bsz, &one, &Hg_sp[bsz2*pij], &bsz, &x_sp[bsz*rj],
                          &inc, &alpha, &v_sp[bsz*ri], &inc);
                    alpha = 1.0;
                }
            }
        }
        break;
    case 1:
//         #pragma omp parallel num_threads(sys.noThreads) private(nri, pii, rj, pij, i, j, nriInt, nb)
        {
//             int this_thread = omp_get_thread_num();
            int this_thread = 0;
//             #pragma omp for
            for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
                nri = ent2entStart[ri+1] - ent2entStart[ri];
                pii = ent2entStart[ri];
                pij = ent2entStart[ri] - 1;

                // Create temporary dense vector:
                for (i = 0; i<nri; i++) {
                    pij += 1;
                    rj = ent2ent[pij];
                    if (rj >= sys.BJ_nrows)
                        break;
                    for (j = 0; j < bsz; j++)
                        sys.xDenseRow_sp[this_thread][bsz*i+j] = x_sp[bsz*rj+j];
                        //SCOPY(&bsz, &x_sp[bsz*rj], &inc, &sys.xDenseRow_sp[this_thread][bsz*i], &inc);
                }

                nriInt = i;
                if (nriInt == 0)
                    error("Error KL6H8P in gmresmpi.cpp was detected.\n");
                
                // Perform matrix-vector product
                nb = bsz*nriInt;
                SGEMV(&chn, &bsz, &nb, &one, &Hg_sp[bsz2*pii], &bsz, &sys.xDenseRow_sp[this_thread][0],
                      &inc, &zero, &v_sp[bsz*ri], &inc);
            }
        }
        break;
    case 2:
        if (bsz != 5)
            error("sys.matvecImplementation == 2 only valid for bsz == 5.\n");
        
        float *Hgpos = &Hg_sp[0];
        float *vpos = &v_sp[0];
        Int *ent2ent_pos = &ent2ent[0];
        
//         Int cacheLineSz = CacheLineSize();
//         Int cacheL0Sz = CacheL0Size();
//         Int cacheL1Sz = CacheL1Size();
//         Int cacheL2Sz = CacheL2Size();
//         Int cacheL3Sz = CacheL3Size();
        
        for (int ri = 0; ri < sys.entpartpts[0]; ri++) {
            int nri = ent2entStart[ri+1] - ent2entStart[ri];
            
            float vtemp1 = 0;
            float vtemp2 = 0;
            float vtemp3 = 0;
            float vtemp4 = 0;
            float vtemp5 = 0;
            
            for (int i=0; i<nri; i++) {
//                 __builtin_prefetch(&x_sp[5 * (*ent2ent_pos + preftchDist)    ] , 0 , 3);
//                 __builtin_prefetch(Hgpos+25*preftchDist      , 0 , 2);
//                 __builtin_prefetch(Hgpos+25*preftchDist + 16 , 0 , 2);
                
                float *xpos = &x_sp[5 * (*ent2ent_pos++)];
                
                vtemp1 += Hgpos[0] * xpos[0] + Hgpos[5] * xpos[1] + Hgpos[10] * xpos[2] + Hgpos[15] * xpos[3] + Hgpos[20] * xpos[4];
                vtemp2 += Hgpos[1] * xpos[0] + Hgpos[6] * xpos[1] + Hgpos[11] * xpos[2] + Hgpos[16] * xpos[3] + Hgpos[21] * xpos[4];
                vtemp3 += Hgpos[2] * xpos[0] + Hgpos[7] * xpos[1] + Hgpos[12] * xpos[2] + Hgpos[17] * xpos[3] + Hgpos[22] * xpos[4];
                vtemp4 += Hgpos[3] * xpos[0] + Hgpos[8] * xpos[1] + Hgpos[13] * xpos[2] + Hgpos[18] * xpos[3] + Hgpos[23] * xpos[4];
                vtemp5 += Hgpos[4] * xpos[0] + Hgpos[9] * xpos[1] + Hgpos[14] * xpos[2] + Hgpos[19] * xpos[3] + Hgpos[24] * xpos[4];
                Hgpos += 25;
                
//                 vtemp1 += (*Hgpos++) * (*xpos);
//                 vtemp2 += (*Hgpos++) * (*xpos);
//                 vtemp3 += (*Hgpos++) * (*xpos);
//                 vtemp4 += (*Hgpos++) * (*xpos);
//                 vtemp5 += (*Hgpos++) * (*xpos);
// 
//                 xpos++;
//                 vtemp1 += (*Hgpos++) * (*xpos);
//                 vtemp2 += (*Hgpos++) * (*xpos);
//                 vtemp3 += (*Hgpos++) * (*xpos);
//                 vtemp4 += (*Hgpos++) * (*xpos);
//                 vtemp5 += (*Hgpos++) * (*xpos);
// 
//                 xpos++;
//                 vtemp1 += (*Hgpos++) * (*xpos);
//                 vtemp2 += (*Hgpos++) * (*xpos);
//                 vtemp3 += (*Hgpos++) * (*xpos);
//                 vtemp4 += (*Hgpos++) * (*xpos);
//                 vtemp5 += (*Hgpos++) * (*xpos);
// 
//                 xpos++;
//                 vtemp1 += (*Hgpos++) * (*xpos);
//                 vtemp2 += (*Hgpos++) * (*xpos);
//                 vtemp3 += (*Hgpos++) * (*xpos);
//                 vtemp4 += (*Hgpos++) * (*xpos);
//                 vtemp5 += (*Hgpos++) * (*xpos);
//                 
//                 xpos++;
//                 vtemp1 += (*Hgpos++) * (*xpos);
//                 vtemp2 += (*Hgpos++) * (*xpos);
//                 vtemp3 += (*Hgpos++) * (*xpos);
//                 vtemp4 += (*Hgpos++) * (*xpos);
//                 vtemp5 += (*Hgpos++) * (*xpos);
            }
            *vpos = vtemp1;
            vpos++;
            *vpos = vtemp2;
            vpos++;
            *vpos = vtemp3;
            vpos++;
            *vpos = vtemp4;
            vpos++;
            *vpos = vtemp5;
            vpos++;
        }
        for (int ri = sys.entpartpts[0]; ri < sys.BJ_nrows; ri++) {
            int nri = ent2entStart[ri+1] - ent2entStart[ri];
            
            float vtemp1 = 0;
            float vtemp2 = 0;
            float vtemp3 = 0;
            float vtemp4 = 0;
            float vtemp5 = 0;
            
            for (int i=0; i<nri; i++) {
                if ((*ent2ent_pos) < sys.BJ_nrows) {
//                 __builtin_prefetch(&x_sp[5 * (*ent2ent_pos + preftchDist)    ] , 0 , 3);
//                 __builtin_prefetch(Hgpos+25*preftchDist      , 0 , 2);
//                 __builtin_prefetch(Hgpos+25*preftchDist + 16 , 0 , 2);
                    
                    float *xpos = &x_sp[5 * (*ent2ent_pos++)];
                    
                    vtemp1 += Hgpos[0] * xpos[0] + Hgpos[5] * xpos[1] + Hgpos[10] * xpos[2] + Hgpos[15] * xpos[3] + Hgpos[20] * xpos[4];
                    vtemp2 += Hgpos[1] * xpos[0] + Hgpos[6] * xpos[1] + Hgpos[11] * xpos[2] + Hgpos[16] * xpos[3] + Hgpos[21] * xpos[4];
                    vtemp3 += Hgpos[2] * xpos[0] + Hgpos[7] * xpos[1] + Hgpos[12] * xpos[2] + Hgpos[17] * xpos[3] + Hgpos[22] * xpos[4];
                    vtemp4 += Hgpos[3] * xpos[0] + Hgpos[8] * xpos[1] + Hgpos[13] * xpos[2] + Hgpos[18] * xpos[3] + Hgpos[23] * xpos[4];
                    vtemp5 += Hgpos[4] * xpos[0] + Hgpos[9] * xpos[1] + Hgpos[14] * xpos[2] + Hgpos[19] * xpos[3] + Hgpos[24] * xpos[4];
                    Hgpos += 25;
                    
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
                }
                else {
                    ent2ent_pos++;
                    Hgpos+=25;
                }
            }
            *vpos = vtemp1;
            vpos++;
            *vpos = vtemp2;
            vpos++;
            *vpos = vtemp3;
            vpos++;
            *vpos = vtemp4;
            vpos++;
            *vpos = vtemp5;
            vpos++;
        }
        break;
    }
    matvectimes[0] += clock() - t;

    /* Wait until all send and receive operations are completely done */
    t = clock();
    MPI_Waitall(*requestCounter, sys.requests, sys.statuses);
    matvectimes[1] += clock() - t;

    /* Copy buffrecv to x_sp */
    for (i=0; i<sys.nentrecv; i++) {
        ri = sys.entrecv[i];
        for (j=0; j<bsz; j++) {
            x_sp[bsz*ri+j] = (float) sys.buffrecv[bsz*i+j];
        }
    }

    /* v <- v + Hg*x for entities NOT in the processor */
    t = clock();
    switch (sys.matvecImplementation) {
    case 0:
//         #pragma omp parallel for num_threads(sys.noThreads) private(nri, rj, pij, i)
        for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
            nri = ent2entStart[ri+1] - ent2entStart[ri];
            for (i = 0; i<nri; i++) {
                pij = ent2entStart[ri]+i;
                rj = ent2ent[pij];
                if (rj >= sys.BJ_nrows) {
                    SGEMV(&chn, &bsz, &bsz, &one, &Hg_sp[bsz2*pij], &bsz, &x_sp[bsz*rj],
                          &inc, &one, &v_sp[bsz*ri], &inc);
                }
            }
        }
        break;
    case 1:
//         #pragma omp parallel num_threads(sys.noThreads) private(nri, pii, rj, pij, pim, i, j, nriExt, nb)
        {
//             int this_thread = omp_get_thread_num();
            int this_thread = 0;
//             #pragma omp for
            for (Int ri = 0; ri < sys.BJ_nrows; ri++) {

                nri = ent2entStart[ri+1] - ent2entStart[ri];
                pij = ent2entStart[ri] - 1;

                for (i = 0; i < nri; i++) {
                    pij += 1;
                    if (ent2ent[pij] >= sys.BJ_nrows)
                        break;
                }
                nriExt = nri - i;

                if (nriExt > 0) {
                    pim = ent2entStart[ri] + nri - nriExt;
                    pij -= 1;

                    // Create temporary dense vector:
                    for (i = 0; i<nriExt; i++) {
                        pij += 1;
                        rj = ent2ent[pij];
                        for (j = 0; j < bsz; j++)
                            sys.xDenseRow_sp[this_thread][bsz*i+j] = x_sp[bsz*rj+j];
                            //SCOPY(&bsz, &x_sp[bsz*rj], &inc, &sys.xDenseRow_sp[this_thread][bsz*i], &inc);
                    }

                    // Perform matrix-vector product:
                    nb = bsz*nriExt;
                    SGEMV(&chn, &bsz, &nb, &one, &Hg_sp[bsz2*pim], &bsz, &sys.xDenseRow_sp[this_thread][0],
                          &inc, &one, &v_sp[bsz*ri], &inc);
                }
            }
        }
        break;
    case 2:
        if (bsz != 5)
            error("sys.matvecImplementation == 2 only valid for bsz == 5.\n");
        
        float *Hgpos = &Hg_sp[25*ent2entStart[sys.entpartpts[0]]];
        float *vpos = &v_sp[5*sys.entpartpts[0]];
        Int *ent2ent_pos = &ent2ent[ent2entStart[sys.entpartpts[0]]];
        
//         Int cacheLineSz = CacheLineSize();
//         Int cacheL0Sz = CacheL0Size();
//         Int cacheL1Sz = CacheL1Size();
//         Int cacheL2Sz = CacheL2Size();
//         Int cacheL3Sz = CacheL3Size();
        
        for (int ri = sys.entpartpts[0]; ri < sys.BJ_nrows; ri++) {
            int nri = ent2entStart[ri+1] - ent2entStart[ri];
            
            float vtemp1 = 0;
            float vtemp2 = 0;
            float vtemp3 = 0;
            float vtemp4 = 0;
            float vtemp5 = 0;
            
            for(int i=0; i<nri; i++) {
                if ((*ent2ent_pos) >= sys.BJ_nrows) {
//                 __builtin_prefetch(&x_sp[5 * (*ent2ent_pos + preftchDist)    ] , 0 , 3);
//                 __builtin_prefetch(Hgpos+25*preftchDist      , 0 , 2);
//                 __builtin_prefetch(Hgpos+25*preftchDist + 16 , 0 , 2);
                    
                    float *xpos = &x_sp[5 * (*ent2ent_pos++)];
                    
                    vtemp1 += Hgpos[0] * xpos[0] + Hgpos[5] * xpos[1] + Hgpos[10] * xpos[2] + Hgpos[15] * xpos[3] + Hgpos[20] * xpos[4];
                    vtemp2 += Hgpos[1] * xpos[0] + Hgpos[6] * xpos[1] + Hgpos[11] * xpos[2] + Hgpos[16] * xpos[3] + Hgpos[21] * xpos[4];
                    vtemp3 += Hgpos[2] * xpos[0] + Hgpos[7] * xpos[1] + Hgpos[12] * xpos[2] + Hgpos[17] * xpos[3] + Hgpos[22] * xpos[4];
                    vtemp4 += Hgpos[3] * xpos[0] + Hgpos[8] * xpos[1] + Hgpos[13] * xpos[2] + Hgpos[18] * xpos[3] + Hgpos[23] * xpos[4];
                    vtemp5 += Hgpos[4] * xpos[0] + Hgpos[9] * xpos[1] + Hgpos[14] * xpos[2] + Hgpos[19] * xpos[3] + Hgpos[24] * xpos[4];
                    Hgpos += 25;
                    
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
                }
                else {
                    ent2ent_pos++;
                    Hgpos+=25;
                }
            }
            *vpos += vtemp1;
            vpos++;
            *vpos += vtemp2;
            vpos++;
            *vpos += vtemp3;
            vpos++;
            *vpos += vtemp4;
            vpos++;
            *vpos += vtemp5;
            vpos++;
        }
        break;
    }
    matvectimes[2] += clock() - t;
    
    /* Convert v to double precision */
    for (i = 0; i < sys.v_sp.size(); i++)
        v[i] = (double) v_sp[i];

    tTotal = clock() - t0;
}

static void matvecMPI_RAS_dp(double* v, sysstruct &sys, double* x, Int * requestCounter, double* matvectimes)
{
    // Description: Double-precision, parallel matrix-vector product 
    //              for RAS(1) preconditioner
    
    Int inc = 1, i, j, n, nri, nriInt, nriExt, pij, pii, pim, ri, rj, rij, nb;
    Int bsz = sys.blkSize;
    Int bsz2    = bsz*bsz;
    Int *ent2entStart = &sys.ent2entStart[0];
    Int *ent2ent = &sys.ent2ent[0];
    double *Hg = &sys.Hg[0];
    double matvectime1, matvectime2, MPIwaitTime, tTotal;
    double one = 1.0, zero = 0.0, alpha;
    char chn = 'N';
    clock_t t, t0 = clock();
    
    /* v <- Hg*x for entities in the processor while sending and receiving data */
    t = clock();
    
    switch (sys.matvecImplementation) {
    case 0:
//         #pragma omp parallel for num_threads(sys.noThreads) private(alpha, nri, rj, pij, i)
        for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
            alpha = 0.0;
            nri = ent2entStart[ri+1] - ent2entStart[ri];
            for (i = 0; i<nri; i++) {
                pij = ent2entStart[ri]+i;
                rj = ent2ent[pij];
                if (rj < sys.BJ_nrows) {
                    DGEMV(&chn, &bsz, &bsz, &one, &Hg[bsz2*pij], &bsz, &x[bsz*rj],
                          &inc, &alpha, &v[bsz*ri], &inc);
                    alpha = 1.0;
                }
            }
        }
        break;
    case 1:
//         #pragma omp parallel num_threads(sys.noThreads) private(nri, pii, rj, pij, i, j, nriInt, nb)
        {
//             int this_thread = omp_get_thread_num();
            int this_thread = 0;
//             #pragma omp for
            for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
                nri = ent2entStart[ri+1] - ent2entStart[ri];
                pii = ent2entStart[ri];
                pij = ent2entStart[ri] - 1;
                
                // Create temporary dense vector:
                for (i = 0; i<nri; i++) {
                    pij += 1;
                    rj = ent2ent[pij];
                    if (rj >= sys.BJ_nrows)
                        break;
                    for (j = 0; j < bsz; j++)
                        sys.xDenseRow[this_thread][bsz*i+j] = x[bsz*rj+j];
                        //DCOPY(&bsz, &x[bsz*rj], &inc, &sys.xDenseRow[this_thread][bsz*i], &inc);
                }

                nriInt = i;
                if (nriInt == 0) {
                    printf("ri = %d, %d, %d, i = %d, nriInt = %d, rj = %d sys.BJ_nrows = %d, pij = %d\n", ri, ent2entStart[ri], ent2entStart[ri+1], i, nriInt, rj, sys.BJ_nrows, pij);
                    error("Error yeahhhh 4FT6H8 in gmresmpi.cpp was detected.\n");
                }
                    
                // Perform matrix-vector product:
                nb = bsz*nriInt;
                DGEMV(&chn, &bsz, &nb, &one, &Hg[bsz2*pii], &bsz, &sys.xDenseRow[this_thread][0],
                      &inc, &zero, &v[bsz*ri], &inc);
            }
        }
        break;
    case 2:
        if (bsz != 5)
            error("sys.matvecImplementation == 2 only valid for bsz == 5.\n");
        
        double *Hgpos = &Hg[0];
        double *vpos = &v[0];
        Int *ent2ent_pos = &ent2ent[0];
        
//         Int cacheLineSz = CacheLineSize();
//         Int cacheL0Sz = CacheL0Size();
//         Int cacheL1Sz = CacheL1Size();
//         Int cacheL2Sz = CacheL2Size();
//         Int cacheL3Sz = CacheL3Size();
        
        for (int ri = 0; ri < sys.entpartpts[0]; ri++) {
            int nri = ent2entStart[ri+1] - ent2entStart[ri];
            
            double vtemp1 = 0;
            double vtemp2 = 0;
            double vtemp3 = 0;
            double vtemp4 = 0;
            double vtemp5 = 0;
            
            for (int i=0; i<nri; i++) {
//                 __builtin_prefetch(&x_sp[5 * (*ent2ent_pos + preftchDist)    ] , 0 , 3);
//                 __builtin_prefetch(Hgpos+25*preftchDist      , 0 , 2);
//                 __builtin_prefetch(Hgpos+25*preftchDist + 16 , 0 , 2);
                
                double *xpos = &x[5 * (*ent2ent_pos++)];
                
                vtemp1 += Hgpos[0] * xpos[0] + Hgpos[5] * xpos[1] + Hgpos[10] * xpos[2] + Hgpos[15] * xpos[3] + Hgpos[20] * xpos[4];
                vtemp2 += Hgpos[1] * xpos[0] + Hgpos[6] * xpos[1] + Hgpos[11] * xpos[2] + Hgpos[16] * xpos[3] + Hgpos[21] * xpos[4];
                vtemp3 += Hgpos[2] * xpos[0] + Hgpos[7] * xpos[1] + Hgpos[12] * xpos[2] + Hgpos[17] * xpos[3] + Hgpos[22] * xpos[4];
                vtemp4 += Hgpos[3] * xpos[0] + Hgpos[8] * xpos[1] + Hgpos[13] * xpos[2] + Hgpos[18] * xpos[3] + Hgpos[23] * xpos[4];
                vtemp5 += Hgpos[4] * xpos[0] + Hgpos[9] * xpos[1] + Hgpos[14] * xpos[2] + Hgpos[19] * xpos[3] + Hgpos[24] * xpos[4];
                Hgpos += 25;
                
//                 vtemp1 += (*Hgpos++) * (*xpos);
//                 vtemp2 += (*Hgpos++) * (*xpos);
//                 vtemp3 += (*Hgpos++) * (*xpos);
//                 vtemp4 += (*Hgpos++) * (*xpos);
//                 vtemp5 += (*Hgpos++) * (*xpos);
// 
//                 xpos++;
//                 vtemp1 += (*Hgpos++) * (*xpos);
//                 vtemp2 += (*Hgpos++) * (*xpos);
//                 vtemp3 += (*Hgpos++) * (*xpos);
//                 vtemp4 += (*Hgpos++) * (*xpos);
//                 vtemp5 += (*Hgpos++) * (*xpos);
// 
//                 xpos++;
//                 vtemp1 += (*Hgpos++) * (*xpos);
//                 vtemp2 += (*Hgpos++) * (*xpos);
//                 vtemp3 += (*Hgpos++) * (*xpos);
//                 vtemp4 += (*Hgpos++) * (*xpos);
//                 vtemp5 += (*Hgpos++) * (*xpos);
// 
//                 xpos++;
//                 vtemp1 += (*Hgpos++) * (*xpos);
//                 vtemp2 += (*Hgpos++) * (*xpos);
//                 vtemp3 += (*Hgpos++) * (*xpos);
//                 vtemp4 += (*Hgpos++) * (*xpos);
//                 vtemp5 += (*Hgpos++) * (*xpos);
//                 
//                 xpos++;
//                 vtemp1 += (*Hgpos++) * (*xpos);
//                 vtemp2 += (*Hgpos++) * (*xpos);
//                 vtemp3 += (*Hgpos++) * (*xpos);
//                 vtemp4 += (*Hgpos++) * (*xpos);
//                 vtemp5 += (*Hgpos++) * (*xpos);
            }
            *vpos = vtemp1;
            vpos++;
            *vpos = vtemp2;
            vpos++;
            *vpos = vtemp3;
            vpos++;
            *vpos = vtemp4;
            vpos++;
            *vpos = vtemp5;
            vpos++;
        }
        for (int ri = sys.entpartpts[0]; ri < sys.BJ_nrows; ri++) {
            int nri = ent2entStart[ri+1] - ent2entStart[ri];
            
            double vtemp1 = 0;
            double vtemp2 = 0;
            double vtemp3 = 0;
            double vtemp4 = 0;
            double vtemp5 = 0;
            
            for (int i=0; i<nri; i++) {
                if ((*ent2ent_pos) < sys.BJ_nrows) {
//                 __builtin_prefetch(&x_sp[5 * (*ent2ent_pos + preftchDist)    ] , 0 , 3);
//                 __builtin_prefetch(Hgpos+25*preftchDist      , 0 , 2);
//                 __builtin_prefetch(Hgpos+25*preftchDist + 16 , 0 , 2);
                    
                    double *xpos = &x[5 * (*ent2ent_pos++)];
                    
                    vtemp1 += Hgpos[0] * xpos[0] + Hgpos[5] * xpos[1] + Hgpos[10] * xpos[2] + Hgpos[15] * xpos[3] + Hgpos[20] * xpos[4];
                    vtemp2 += Hgpos[1] * xpos[0] + Hgpos[6] * xpos[1] + Hgpos[11] * xpos[2] + Hgpos[16] * xpos[3] + Hgpos[21] * xpos[4];
                    vtemp3 += Hgpos[2] * xpos[0] + Hgpos[7] * xpos[1] + Hgpos[12] * xpos[2] + Hgpos[17] * xpos[3] + Hgpos[22] * xpos[4];
                    vtemp4 += Hgpos[3] * xpos[0] + Hgpos[8] * xpos[1] + Hgpos[13] * xpos[2] + Hgpos[18] * xpos[3] + Hgpos[23] * xpos[4];
                    vtemp5 += Hgpos[4] * xpos[0] + Hgpos[9] * xpos[1] + Hgpos[14] * xpos[2] + Hgpos[19] * xpos[3] + Hgpos[24] * xpos[4];
                    Hgpos += 25;
                    
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
                }
                else {
                    ent2ent_pos++;
                    Hgpos+=25;
                }
            }
            *vpos = vtemp1;
            vpos++;
            *vpos = vtemp2;
            vpos++;
            *vpos = vtemp3;
            vpos++;
            *vpos = vtemp4;
            vpos++;
            *vpos = vtemp5;
            vpos++;
        }
        break;
    }
    matvectimes[0] += clock() - t;
    
    /* Wait until all send and receive operations are completely done */
    t = clock();
    MPI_Waitall(*requestCounter, sys.requests, sys.statuses);
    matvectimes[1] += clock() - t;
    
    /* Copy buffrecv to x */
    for (i=0; i<sys.nentrecv; i++) {
        ri = sys.entrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &x[bsz*ri], &inc);
    }

    /* v <- v + Hg*x for entities NOT in the processor */
    t = clock();
    switch (sys.matvecImplementation) {
    case 0:
//         #pragma omp parallel for num_threads(sys.noThreads) private(nri, rj, pij, i)
        for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
            nri = ent2entStart[ri+1] - ent2entStart[ri];
            for (i = 0; i<nri; i++) {
                pij = ent2entStart[ri]+i;
                rj = ent2ent[pij];
                if (rj >= sys.BJ_nrows) {
                    DGEMV(&chn, &bsz, &bsz, &one, &Hg[bsz2*pij], &bsz, &x[bsz*rj],
                          &inc, &one, &v[bsz*ri], &inc);
                }
            }
        }
        break;
    case 1:
//         #pragma omp parallel num_threads(sys.noThreads) private(nri, pii, rj, pij, pim, i, j, nriExt, nb)
        {
//             int this_thread = omp_get_thread_num();
            int this_thread = 0;
//             #pragma omp for
            for (Int ri = 0; ri < sys.BJ_nrows; ri++) {
                nri = ent2entStart[ri+1] - ent2entStart[ri];
                pij = ent2entStart[ri] - 1;

                for (i = 0; i < nri; i++) {
                    pij += 1;
                    if (ent2ent[pij] >= sys.BJ_nrows)
                        break;
                }
                nriExt = nri - i;

                if (nriExt > 0) {
                    pim = ent2entStart[ri] + nri - nriExt;
                    pij -= 1;

                    // Create temporary dense vector:
                    for (i = 0; i<nriExt; i++) {
                        pij += 1;
                        rj = ent2ent[pij];
                        for (j = 0; j < bsz; j++)
                            sys.xDenseRow[this_thread][bsz*i+j] = x[bsz*rj+j];
                            //DCOPY(&bsz, &x[bsz*rj], &inc, &xDenseRow[this_thread][bsz*i], &inc);
                    }

                    // Perform matrix-vector product:
                    nb = bsz*nriExt;
                    DGEMV(&chn, &bsz, &nb, &one, &Hg[bsz2*pim], &bsz, &sys.xDenseRow[this_thread][0],
                          &inc, &one, &v[bsz*ri], &inc);
                }
            }
        }
        break;
    case 2:
        if (bsz != 5)
            error("sys.matvecImplementation == 2 only valid for bsz == 5.\n");
        
        double *Hgpos = &Hg[25*ent2entStart[sys.entpartpts[0]]];
        double *vpos = &v[5*sys.entpartpts[0]];
        Int *ent2ent_pos = &ent2ent[ent2entStart[sys.entpartpts[0]]];
        
//         Int cacheLineSz = CacheLineSize();
//         Int cacheL0Sz = CacheL0Size();
//         Int cacheL1Sz = CacheL1Size();
//         Int cacheL2Sz = CacheL2Size();
//         Int cacheL3Sz = CacheL3Size();
        
        for (int ri = sys.entpartpts[0]; ri < sys.BJ_nrows; ri++) {
            int nri = ent2entStart[ri+1] - ent2entStart[ri];
            
            double vtemp1 = 0;
            double vtemp2 = 0;
            double vtemp3 = 0;
            double vtemp4 = 0;
            double vtemp5 = 0;
            
            for(int i=0; i<nri; i++) {
                if ((*ent2ent_pos) >= sys.BJ_nrows) {
//                 __builtin_prefetch(&x_sp[5 * (*ent2ent_pos + preftchDist)    ] , 0 , 3);
//                 __builtin_prefetch(Hgpos+25*preftchDist      , 0 , 2);
//                 __builtin_prefetch(Hgpos+25*preftchDist + 16 , 0 , 2);
                    
                    double *xpos = &x[5 * (*ent2ent_pos++)];
                    
                    vtemp1 += Hgpos[0] * xpos[0] + Hgpos[5] * xpos[1] + Hgpos[10] * xpos[2] + Hgpos[15] * xpos[3] + Hgpos[20] * xpos[4];
                    vtemp2 += Hgpos[1] * xpos[0] + Hgpos[6] * xpos[1] + Hgpos[11] * xpos[2] + Hgpos[16] * xpos[3] + Hgpos[21] * xpos[4];
                    vtemp3 += Hgpos[2] * xpos[0] + Hgpos[7] * xpos[1] + Hgpos[12] * xpos[2] + Hgpos[17] * xpos[3] + Hgpos[22] * xpos[4];
                    vtemp4 += Hgpos[3] * xpos[0] + Hgpos[8] * xpos[1] + Hgpos[13] * xpos[2] + Hgpos[18] * xpos[3] + Hgpos[23] * xpos[4];
                    vtemp5 += Hgpos[4] * xpos[0] + Hgpos[9] * xpos[1] + Hgpos[14] * xpos[2] + Hgpos[19] * xpos[3] + Hgpos[24] * xpos[4];
                    Hgpos += 25;
                    
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
// 
//                     xpos++;
//                     vtemp1 += (*Hgpos++) * (*xpos);
//                     vtemp2 += (*Hgpos++) * (*xpos);
//                     vtemp3 += (*Hgpos++) * (*xpos);
//                     vtemp4 += (*Hgpos++) * (*xpos);
//                     vtemp5 += (*Hgpos++) * (*xpos);
                }
                else {
                    ent2ent_pos++;
                    Hgpos+=25;
                }
            }
            *vpos += vtemp1;
            vpos++;
            *vpos += vtemp2;
            vpos++;
            *vpos += vtemp3;
            vpos++;
            *vpos += vtemp4;
            vpos++;
            *vpos += vtemp5;
            vpos++;
        }
        break;
    }
    matvectimes[2] += clock() - t;

    tTotal = clock() - t0;
}

static void matvecMPI_RAS(double* v, sysstruct &sys, double* x, Int * requestCounter, double* matvectimes)
{
    if (sys.matvecPrecision == 0 && sys.robustMode == 0)
        matvecMPI_RAS_sp(v, sys, x, requestCounter, &matvectimes[0]);
    else
        matvecMPI_RAS_dp(v, sys, x, requestCounter, &matvectimes[0]);
}

void matvecMPI(double* v, sysstruct &sys, double* x, Int sendInputVector, Int sendOutputVector, Int * requestCounter, double* matvectimes)
{
    // Send input vector (if necessary). Reception is done in matvec
    if (sendInputVector == 1 && sys.matvecImplementation != 2)
        sendvector(sys, x, requestCounter);
    
    // Compute matrix-vector product:
    matvecMPI_RAS(v, sys, x, requestCounter, &matvectimes[0]);
    
    // Send output vector (if necessary). Reception will be done in other part of the code
    if (sendOutputVector == 1)
        sendvector(sys, v, requestCounter);
}

#endif

//static void matvec_v1(double* v, sysstruct &sys, double* x)
//{
//    char chn = 'N';
//    double one = 1.0, zero = 0.0;
//
//    Int inc = 1, i, nri, pij = -1, rj;
//
//    Int nrows   = sys.numEntities;
//    Int *rowpts = &sys.ent2entStart[0];
//    Int *colind = &sys.ent2ent[0];
//    Int bsz     = sys.blkSize;
//    Int bsz2    = bsz*bsz;
//
//    for (Int ri = 0; ri < nrows; ri++) {
//        nri = rowpts[ri+1] - rowpts[ri];
//        pij = rowpts[ri];
//        rj = colind[pij];
//        DGEMV(&chn, &bsz, &bsz, &one, &sys.Hg[bsz2*pij], &bsz, &x[bsz*rj],
//                  &inc, &zero, &v[bsz*ri], &inc);
//        for (i = 1; i<nri; i++) {
//            pij += 1;
//            rj = colind[pij];
//            DGEMV(&chn, &bsz, &bsz, &one, &sys.Hg[bsz2*pij], &bsz, &x[bsz*rj],
//                  &inc, &one, &v[bsz*ri], &inc);
//        }
//    }
//}

static void matvec_dp(double* v, sysstruct &sys, double* x)
{
    char chn = 'N';
    double one = 1.0, zero = 0.0, alpha;

    Int inc = 1, i, nri, pij, rj;
    Int nrows   = sys.numEntities;
    Int *rowpts = &sys.ent2entStart[0];
    Int *colind = &sys.ent2ent[0];
    Int bsz     = sys.blkSize;
    Int bsz2    = bsz*bsz;

    //pij = -1;
//     #pragma omp parallel for num_threads(sys.noThreads) private(alpha, nri, rj, pij, i)
    for (Int ri = 0; ri < nrows; ri++) {
        alpha = 0.0;
        nri = rowpts[ri+1] - rowpts[ri];
        for (i = 0; i<nri; i++) {
            //pij += 1;
            pij = rowpts[ri]+i;
            rj = colind[pij];
            DGEMV(&chn, &bsz, &bsz, &one, &sys.Hg[bsz2*pij], &bsz, &x[bsz*rj],
                  &inc, &alpha, &v[bsz*ri], &inc);
            alpha = 1.0;
        }
    }
}

static void matvec_sp(double* v, sysstruct &sys, double* x)
{
    char chn = 'N';
    float one = 1.0, zero = 0.0, alpha;
    float *Hg_sp = &sys.Hg_sp[0];
    float *x_sp = &sys.x_sp[0];
    float *v_sp = &sys.v_sp[0];
    
    Int inc = 1, i, nri, pij, rj;
    Int nrows   = sys.numEntities;
    Int *rowpts = &sys.ent2entStart[0];
    Int *colind = &sys.ent2ent[0];
    Int bsz     = sys.blkSize;
    Int bsz2    = bsz*bsz;

    /* Convert x to single precision */
    for (i = 0; i < sys.x_sp.size(); i++)
        x_sp[i] = (float) x[i];
    
    //pij = -1;
//     #pragma omp parallel for num_threads(sys.noThreads) private(alpha, nri, rj, pij, i)
    for (Int ri = 0; ri < nrows; ri++) {
        alpha = 0.0;
        nri = rowpts[ri+1] - rowpts[ri];
        for (i = 0; i<nri; i++) {
            //pij += 1;
            pij = rowpts[ri]+i;
            rj = colind[pij];
            SGEMV(&chn, &bsz, &bsz, &one, &sys.Hg_sp[bsz2*pij], &bsz, &x_sp[bsz*rj],
                  &inc, &alpha, &v_sp[bsz*ri], &inc);
            alpha = 1.0;
        }
    }
    
    /* Convert v to double precision */
    for (i = 0; i < sys.v_sp.size(); i++)
        v[i] = (double) v_sp[i];
}

void matvec(double* v, sysstruct &sys, double* x)
{
    if (sys.matvecPrecision == 0 && sys.robustMode == 0)
        matvec_sp(v, sys, x);
    else
        matvec_dp(v, sys, x);
}

//static void matvec_v3(double* v, sysstruct &sys, double* x)
//{
//    char chn = 'N';
//    double one = 1.0, zero = 0.0;
//
//    Int inc = 1, i, j, nri, pii, pij, rj;
//
//    Int nrows   = sys.numEntities;
//    Int *rowpts = &sys.ent2entStart[0];
//    Int *colind = &sys.ent2ent[0];
//    Int bsz     = sys.blkSize;
//    Int bsz2    = bsz*bsz;
//    Int na = bsz, nb;
//
//    for (Int ri = 0; ri < nrows; ri++) {
//        nri = rowpts[ri+1] - rowpts[ri];
//        pii = rowpts[ri];
//        ri = colind[pii];
//
//        // Create temporary dense vector:
//        for (i = 0; i<nri; i++) {
//            pij = pii+i;
//            rj = colind[pij];
//            for (j = 0; j < bsz; j++)
//                 sys.xDenseRow[bsz*i+j] = x[bsz*rj+j];
//            //DCOPY(&bsz, &x[bsz*rj], &inc, &sys.xDenseRow[bsz*i], &inc);
//        }
//
//        nb = bsz*nri;
//        DGEMV(&chn, &na, &nb, &one, &sys.Hg[bsz2*pii], &na, &sys.xDenseRow[0],
//              &inc, &zero, &v[bsz*ri], &inc);
//    }
//}

//static void matvec_v4(double* v, sysstruct &sys, double* x)
//{
//    char chn = 'N';
//    double one = 1.0, zero = 0.0;
//
//    Int inc = 1, i, j, nri, pii, pij, rj;
//
//    Int nrows   = sys.numEntities;
//    Int *rowpts = &sys.ent2entStart[0];
//    Int *colind = &sys.ent2ent[0];
//    Int bsz     = sys.blkSize;
//    Int bsz2    = bsz*bsz;
//    Int na = bsz, nb;
//
//    // Create dense vector:
//    for (Int ri = 0; ri < nrows; ri++) {
//        nri = rowpts[ri+1] - rowpts[ri];
//        pii = rowpts[ri];
//
//        for (i = 0; i<nri; i++) {
//            pij = pii+i;
//            rj = colind[pij];
//            DCOPY(&bsz, &x[bsz*rj], &inc, &sys.xDense[bsz*pij], &inc);
//        }
//    }
//
//    // Perform matrix-vector product
//    for (Int ri = 0; ri < nrows; ri++) {
//        nri = rowpts[ri+1] - rowpts[ri];
//        pii = rowpts[ri];
//        ri = colind[pii];
//
//        nb = bsz*nri;
//        DGEMV(&chn, &na, &nb, &one, &sys.Hg[bsz2*pii], &na, &sys.xDense[bsz*pii],
//              &inc, &zero, &v[bsz*ri], &inc);
//    }
//}
//
//static void matvec_v5(double* v, sysstruct &sys, double* x)
//{
//    // Note: This implementation showed very slow performance
//    char chn = 'N';
//    double one = 1.0, zero = 0.0;
//
//    Int inc = 1, i, nri, pij = -1, rj, k, kk;
//
//    Int nrows   = sys.numEntities;
//    Int *rowpts = &sys.ent2entStart[0];
//    Int *colind = &sys.ent2ent[0];
//    Int bsz     = sys.blkSize;
//    Int bsz2    = bsz*bsz;
//
//    for (Int ri = 0; ri < nrows; ri++) {
//        nri = rowpts[ri+1] - rowpts[ri];
//        for (k=0; k<bsz; k++) {
//            v[bsz*ri+k] = 0;
//        }
//        for (i = 0; i<nri; i++) {
//            pij += 1;
//            rj = colind[pij];
//            for (kk=0; kk<bsz; kk++) {
//                for (k=0; k<bsz; k++) {
//                    v[bsz*ri+k] += sys.Hg[bsz2*pij+k+bsz*kk]*x[bsz*rj+kk];
//                }
//            }
//        }
//    }
//}

#endif
