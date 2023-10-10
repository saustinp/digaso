#ifndef __PSLAK
#define __PSLAK

// Written by: C. Nguyen & P. Fernandez

#ifdef  HAVE_MPI

vector <Int> sys.ent2OwnEntEndJ[numEntities] // Contains the final location, so that "nriOwn = ent2OwnEntEndJ[ri] - ent2entStartJ[ri] + 1"
//vector <Int> sys.buffrecvStart[sys.nentrecv+1];
        
void matvecfreeMPI(double* v, sysstruct &sys, double* x, Int sendInputVector)
{
    error("(1) Validate and (2) thread matvecfreeMPI.\n");
    
    // Note 1: Same implementation both for single and double precision matrix-vector product
    // Note 2: The current implementation is a bit hacky. udg and uh are used as 
    //     temporary vectors to store sol.UDG and sol.UH since computeSchurResidualMPI
    //     uses sol.UDG and sol.UH to compute the new residual.
    // Note 3: sys.x may have been sent somewhere in the code and may need to be received here.

    Int ne = ndims[5];
    Int npv = ndims[9];
    Int ncu = ndims[20];
    Int n1 = npv*ncu;
    
    double *udg = &sol.UDG4FD[0];
    double *uh = &sol.UH4FD[0];
    double *duh = &x[0];
    double *Rg1 = &sys.Rg4FD[0];
    double *r1 = &sys.r4FD[0];
    double *Rg0 = &sys.Rg_0[0];
    double epsilon = 1.0e-6;
    
    for (Int i = 0; i < sol.UDG.size(); i++)
        udg[i] = sol.UDG[i];
    for (Int i = 0; i < sol.UH.size(); i++)
        uh[i] = sol.UH[i];
    
    if (sendInputVector == 0)
        updateSol(sys, elems, mesh, master, sol, app, temps, ndims, &sol.UDG[0], &sol.UH[0], duh, epsilon);
    else if (sendInputVector == 1)
        updateSolMPI(sys, elems, mesh, master, sol, app, temps, ndims, &sol.UDG[0], &sol.UH[0], duh, epsilon);
    
    computeSchurResidual(sys, elems, mesh, master, sol, app, temps, ndims, &Rg1[0], &r1[0]);
    
    for (Int i = 0; i < sol.UDG.size(); i++)
        sol.UDG[i] = udg[i];
    for (Int i = 0; i < sol.UH.size(); i++)
        sol.UH[i] = uh[i];
    
    Int n1 = sys.rowStartR[sys.numIntEnt];
    for (Int i = 0; i < n1; i++)
        v[i] = (Rg1[i] - Rg0[i]) / epsilon;
}

void matvecMPI_sp(double* v, sysstruct &sys, double* x, Int * requestCounter, double* matvectimes)
{
    // Single-precision parallel matrix-vector product: v = Hg * x
    
    Int inc = 1;
    Int numOwnEnt = sys.numOwnEnt;
    Int *ent2entStartJ = &sys.ent2entStartJ[0];
    Int *ent2entJ = &sys.ent2entJ[0];
    Int *ent2OwnEntEndJ = &sys.ent2OwnEntEndJ[0];
    Int *rowStartJ = &sys.rowStartJ[0];
    Int *blockStartJ = &sys.blockStartJ[0];
    Int *rowStartR = &sys.rowStartR[0];
    Int *entSz = &sys.entSz[0];
    //Int *buffrecvStart = &sys.buffrecvStart[0];
    float *Hg_sp = &sys.Hg_sp[0];
    float *x_sp = &sys.x_sp[0];
    float *v_sp = &sys.v_sp[0];
    double matvectime1, matvectime2, MPIwaitTime, tTotal;
    float one = 1.0, zero = 0.0, alpha;
    char chn = 'N';
    clock_t t, t0 = clock();
    
    // Convert x to single precision:
    sys.x_sp.insert(sys.x_sp.begin(), sys.x.begin(), sys.x.end());
    
    // v <- Hg*x for entities in the processor while sending and receiving data:
    t = clock();
    switch (sys.matvecImplementation) {
    case 0:
//         #pragma omp parallel for num_threads(sys.noThreads) //private(alpha, nri, rj, pij, i)
        for (Int ri = 0; ri < numOwnEnt; ri++) {
            double alpha = 0.0;
            Int bsz_row = entSz[ri];
            Int nriOwn = ent2OwnEntEndJ[ri] - ent2entStartJ[ri] + 1;
            Int pij = ent2entStartJ[ri];
            Int Hg_pts = rowStartJ[ri];
            Int v_pts = rowStartR[ri];
            
            if (nriOwn <= 0)
                error("Error AEF7K9 in matvecMPI_sp was detected.\n");
            
            for (Int i = 0; i < nriOwn; i++) {
                Int rj = ent2entJ[pij];
                Int bsz_col = sys.entSz[rj];
                SGEMV(&chn, &bsz_row, &bsz_col, &one, &Hg_sp[Hg_pts], &bsz_row, &x_sp[rowStartR[rj]],
                      &inc, &alpha, &v_sp[v_pts], &inc);
                pij++;
                Hg_pts += bsz_row*bsz_col;
                alpha = 1.0;
            }
        }
        break;
    case 1:
//         #pragma omp parallel num_threads(sys.noThreads) //private(nri, pii, rj, pij, i, j, nriInt, nb)
        {
//             int this_thread = omp_get_thread_num();
            int this_thread = 0;
//             #pragma omp for
            for (Int ri = 0; ri < numOwnEnt; ri++) {
                Int bsz_row = sys.entSz[ri];
                Int nriOwn = ent2OwnEntEndJ[ri] - ent2entStartJ[ri] + 1;
                Int pij = ent2entStartJ[ri];
                
                if (nriOwn <= 0)
                    error("Error AEF7K9 in matvecMPI_sp was detected.\n");
                
                // Create temporary dense vector:
                std::vector<float>::iterator it = sys.xDenseRow_sp[this_thread].begin();
                Int bsz_cols = 0;
                for (Int i = 0; i < nriOwn; i++) {
                    Int rj = ent2entJ[pij++];
                    bsz_cols += sys.entSz[rj];//sys.entStart[rj+1] - sys.entStart[rj];
                    it = sys.xDenseRow_sp[this_thread].insert(it, x_sp+rowStartR[rj], x_sp+rowStartR[rj+1]);
//                     SCOPY(&bsz, &x_sp[bsz*rj], &inc, &sys.xDenseRow_sp[this_thread][bsz*i], &inc);
                }
                
                // Perform matrix-vector product:
                SGEMV(&chn, &bsz_row, &bsz_cols, &one, &Hg_sp[rowStartJ[ri]], &bsz_row, &sys.xDenseRow_sp[this_thread][0],
                      &inc, &zero, &v_sp[rowStartR[ri]], &inc);
            }
        }
        break;
    }
    matvectimes[0] += clock() - t;
    
    /* Wait until all send and receive operations are completely done */
    t = clock();
    MPI_Waitall(*requestCounter, sys.requests, sys.statuses);
    matvectimes[1] += clock() - t;
    
    /* Copy buffrecv to x_sp */
    Int buff_pt = 0;
    for (Int i=0; i<sys.nentrecv; i++) {
        Int ri = sys.entrecv[i];
        Int x_pt = rowStartR[ri];
        //Int buff_pt = buffrecvStart[i];
        Int len = entSz[ri];
        for (Int j=0; j<len; j++) {
            x_sp[x_pt++] = (float) sys.buffrecv[buff_pt++];
        }
    }
    
    /* v <- v + Hg*x for entities NOT in the processor */
    t = clock();
    switch (sys.matvecImplementation) {
    case 0:
//         #pragma omp parallel for num_threads(sys.noThreads) //private(nri, rj, pij, i)
        for (Int ri = 0; ri < numOwnEnt; ri++) {
            Int bsz_row = entSz[ri];
            Int pij = ent2OwnEntEndJ[ri] + 1;
            Int nriNotOwn = ent2entStartJ[ri+1] - ent2OwnEntEndJ[ri] - 1;
            Int Hg_pts = rowStartJ[ri];
            Int v_pts = rowStartR[ri];
            
            if (nriNotOwn < 0)
                error("Error KL6H8P in matvecMPI_sp was detected.\n");
            
            for (Int i = 0; i < nriNotOwn; i++) {
                Int rj = ent2entJ[pij];
                Int bsz_col = sys.entSz[rj];
                SGEMV(&chn, &bsz_row, &bsz_col, &one, &Hg_sp[Hg_pts], &bsz_row, &x_sp[rowStartR[rj]],
                      &inc, &one, &v_sp[v_pts], &inc);
                pij++;
                Hg_pts += bsz_row*bsz_col;
            }
        }
        break;
    case 1:
//         #pragma omp parallel num_threads(sys.noThreads) //private(nri, pii, rj, pij, pim, i, j, nriExt, nb)
        {
//             int this_thread = omp_get_thread_num();
            int this_thread = 0;
//             #pragma omp for
            for (Int ri = 0; ri < numOwnEnt; ri++) {
                Int bsz_row = sys.entSz[ri];
                Int pij = ent2OwnEntEndJ[ri] + 1;
                Int nriNotOwn = ent2entStartJ[ri+1] - ent2OwnEntEndJ[ri] - 1;
                if (nriNotOwn < 0)
                    error("Error KL6H8P in matvecMPI_sp was detected.\n");
                
                
                // Create temporary dense vector:
                std::vector<float>::iterator it = sys.xDenseRow_sp[this_thread].begin();
                Int bsz_cols = 0;
                for (Int i = 0; i < nriNotOwn; i++) {
                    Int rj = ent2entJ[pij++];
                    bsz_cols += sys.entSz[rj];
                    it = sys.xDenseRow_sp[this_thread].insert(it, x_sp+rowStartR[rj], x_sp+rowStartR[rj+1]);
//                     SCOPY(&bsz, &x_sp[bsz*rj], &inc, &sys.xDenseRow_sp[this_thread][bsz*i], &inc);
                }

                // Perform matrix-vector product:
                SGEMV(&chn, &bsz_row, &bsz_cols, &one, &Hg_sp[rowStartJ[ri]], &bsz_row, &sys.xDenseRow_sp[this_thread][0],
                      &inc, &one, &v_sp[rowStartR[ri]], &inc);
            }
        }
        break;
    }
    matvectimes[2] += clock() - t;
    
    /* Convert v to double precision */
    sys.v.insert(sys.v.begin(), sys.v_sp.begin(), sys.v_sp.end());
    
    tTotal = clock() - t0;
}

void matvecMPI_dp(double* v, sysstruct &sys, double* x, Int * requestCounter, double* matvectimes)
{
    // Double-precision, parallel matrix-vector product
    
    MODIFY FROM _sp ONCE _sp IS VALIDATED
    
    Int i, j, n, pij, pii, pim, ri, rij, inc = 1;
    Int *ent2entStartJ = &sys.ent2entStartJ[0];
    Int *ent2entJ = &sys.ent2entJ[0];
    Int *ent2OwnEntEndJ = &sys.ent2OwnEntEndJ[0];
    Int *blockStartJ = &sys.blockStartJ[0];
    Int *blockStartR = &sys.blockStartR[0];
    Int *entSz = &sys.entSz[0];
    double *Hg = &sys.Hg[0];
    double matvectime1, matvectime2, MPIwaitTime, tTotal;
    double one = 1.0, zero = 0.0, alpha;
    char chn = 'N';
    clock_t t, t0 = clock();
    
    /* v <- Hg*x for entities in the processor while sending and receiving data */
    t = clock();
    switch (sys.matvecImplementation) {
    case 0:
//         #pragma omp parallel for num_threads(sys.noThreads) //private(alpha, nri, rj, pij, i)
        for (Int ri = 0; ri < numOwnEnt; ri++) {
            double alpha = 0.0;
            Int bsz_row = entSz[ri];
            Int nriOwn = ent2OwnEntEndJ[ri] - ent2entStartJ[ri];
            Int pij = ent2entStartJ[ri];
            if (nriOwn <= 0)
                error("Error KL6H8P in gmresmpi.cpp was detected.\n");
            
            for (Int i = 0; i < nriOwn; i++) {
                Int rj = ent2entJ[pij];
                Int bsz_col = sys.entSz[rj];
                DGEMV(&chn, &bsz_row, &bsz_col, &one, &Hg[blockStartJ[pij]], &bsz_row, &x[blockStartR[rj]],
                      &inc, &alpha, &v[blockStartR[ri]], &inc);
                pij++;
                alpha = 1.0;
            }
        }
        break;
    case 1:
//         #pragma omp parallel num_threads(sys.noThreads) //private(nri, pii, rj, pij, i, j, nriInt, nb)
        {
//             int this_thread = omp_get_thread_num();
            int this_thread = 0;
//             #pragma omp for
            for (Int ri = 0; ri < numOwnEnt; ri++) {
                Int bsz_row = sys.entSz[ri];
                Int nriOwn = ent2OwnEntEndJ[ri] - ent2entStartJ[ri];
                Int pij = ent2entStartJ[ri];
                if (nriOwn <= 0)
                    error("Error KL6H8P in gmresmpi.cpp was detected.\n");
                
                // Create temporary dense vector:
                std::vector<double>::iterator it = sys.xDenseRow[this_thread].begin();
                Int bsz_col = 0;
                for (Int i = 0; i < nriOwn; i++) {
                    Int rj = ent2entJ[pij++];
                    bsz_col += sys.entStart[rj+1] - sys.entStart[rj];
                    it = sys.xDenseRow[this_thread].insert(it, x+blockStartR[rj], x+blockStartR[rj+1]);
//                     DCOPY(&bsz, &x[bsz*rj], &inc, &sys.xDenseRow[this_thread][bsz*i], &inc);
                }

                // Perform matrix-vector product:
                DGEMV(&chn, &bsz_row, &bsz_col, &one, &Hg[blockStartJ[ent2entStartJ[ri]]], &bsz_row, &sys.xDenseRow[this_thread][0],
                      &inc, &zero, &v[blockStartR[ri]], &inc);
            }
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
        DCOPY(&entSz[ri], &sys.buffrecv[buffrecvStart[i]], &inc, &x[blockStartR[ri]], &inc);
//         for (j=0; j<entSz[ri]; j++) {
//             x[blockStartR[ri]+j] = sys.buffrecv[buffrecvStart[i]+j];
//         }
    }
    
    /* v <- v + Hg*x for entities NOT in the processor */
    t = clock();
    switch (sys.matvecImplementation) {
    case 0:
//         #pragma omp parallel for num_threads(sys.noThreads) //private(nri, rj, pij, i)
        for (Int ri = 0; ri < numOwnEnt; ri++) {
            Int nriNotOwn = ent2entStartJ[ri+1] - ent2OwnEntEndJ[ri];
            Int pij = ent2OwnEntEndJ[ri];
            Int bsz_row = sys.entSz[ri];
            if (nriNotOwn < 0)
                error("Error KL6H8P in gmresmpi.cpp was detected.\n");
            
            for (Int i = 0; i < nriNotOwn; i++) {
                Int rj = ent2entJ[pij];
                Int bsz_col = sys.entSz[rj];
                DGEMV(&chn, &bsz_row, &bsz_col, &one, &Hg[blockStartJ[pij]], &bsz_row, &x[blockStartR[rj]],
                      &inc, &one, &v[blockStartR[ri]], &inc);
                pij++;
            }
        }
        break;
    case 1:
//         #pragma omp parallel num_threads(sys.noThreads) //private(nri, pii, rj, pij, pim, i, j, nriExt, nb)
        {
//             int this_thread = omp_get_thread_num();
            int this_thread = 0;
//             #pragma omp for
            for (Int ri = 0; ri < numOwnEnt; ri++) {
                Int bsz_row = sys.entSz[ri];
                Int nriNotOwn = ent2entStartJ[ri+1] - ent2OwnEntEndJ[ri];
                Int pij = ent2OwnEntEndJ[ri];
                if (nriNotOwn < 0)
                    error("Error KL6H8P in gmresmpi.cpp was detected.\n");
                
                // Create temporary dense vector:
                std::vector<double>::iterator it = sys.xDenseRow[this_thread].begin();
                Int bsz_col = 0;
                for (Int i = 0; i < nriNotOwn; i++) {
                    Int rj = ent2entJ[pij++];
                    bsz_col += sys.entStart[rj+1] - sys.entStart[rj];
                    it = sys.xDenseRow[this_thread].insert(it, x+blockStartR[rj], x+blockStartR[rj+1]);
//                     DCOPY(&bsz, &x[bsz*rj], &inc, &sys.xDenseRow[this_thread][bsz*i], &inc);
                }

                // Perform matrix-vector product:
                DGEMV(&chn, &bsz_row, &bsz_col, &one, &Hg[blockStartJ[ent2OwnEntEndJ[ri]]], &bsz_row, &sys.xDenseRow[this_thread][0],
                      &inc, &one, &v[blockStartR[ri]], &inc);
            }
        }
        break;
    }
    matvectimes[2] += clock() - t;
    
    tTotal = clock() - t0;
}

void matvecMPI(double* v, sysstruct &sys, double* x, Int sendInputVector, Int sendOutputVector, Int * requestCounter, double* matvectimes)
{
    // Send input vector (if necessary). Reception is done in matvec
    if (sendInputVector == 1 && ?? sys.matvecImplementation != 2)
        sendvector(sys, x, requestCounter);
    
    // Compute matrix-vector product:
    if (sys.matvecImplementation == -1)
        matvecfreeMPI(v, sys, x, sendInputVector);
    else if (sys.matvecPrecision == 0 && sys.robustMode == 0)
        matvecMPI_sp(v, sys, x, requestCounter, &matvectimes[0]);
    else
        matvecMPI_dp(v, sys, x, requestCounter, &matvectimes[0]);

    // Send output vector (if necessary). Reception will be done in other part of the code
    if (sendOutputVector == 1)
        sendvector(sys, v, requestCounter);
}

#endif

//void matvec_v1(double* v, sysstruct &sys, double* x)
//{
//    char chn = 'N';
//    double one = 1.0, zero = 0.0;
//
//    Int inc = 1, i, nri, pij = -1, rj;
//
//    Int nrows   = sys.numEntities;
//    Int *rowpts = &sys.ent2entStartJ[0];
//    Int *colind = &sys.ent2entJ[0];
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

void matvec_dp(double* v, sysstruct &sys, double* x)
{
    error("matvec_dp not developed yet for _v1.0 of the code.\n");
    
    char chn = 'N';
    double one = 1.0, zero = 0.0, alpha;
    
    Int inc = 1, i, nri, pij, rj;
    Int nrowsJ   = sys.numEntitiesJ;
    Int *rowptsJ = &sys.ent2entStartJ[0];
    Int *colindJ = &sys.ent2entJ[0];
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

void matvec_sp(double* v, sysstruct &sys, double* x)
{
    error("matvec_sp not developed yet for _v1.0 of the code.\n");
    
    char chn = 'N';
    float one = 1.0, zero = 0.0, alpha;
    float *Hg_sp = &sys.Hg_sp[0];
    float *x_sp = &sys.x_sp[0];
    float *v_sp = &sys.v_sp[0];
    
    Int inc = 1, i, nri, pij, rj;
    Int nrowsJ   = sys.numEntitiesJ;
    Int *rowptsJ = &sys.ent2entStartJ[0];
    Int *colindJ = &sys.ent2entJ[0];
    Int bsz     = sys.blkSize;
    Int bsz2    = bsz*bsz;

    /* Convert x to single precision */
//     for (i = 0; i < sys.x_sp.size(); i++)
//         x_sp[i] = (float) x[i];
    sys.x_sp.insert(sys.x_sp.begin(), sys.x.begin(), sys.x.end());
    
    //pij = -1;
//     #pragma omp parallel for num_threads(sys.noThreads) private(alpha, nri, rj, pij, i)
    for (Int ri = 0; ri < nrowsJ; ri++) {
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
//     for (i = 0; i < sys.v_sp.size(); i++)
//         v[i] = (double) v_sp[i];
    sys.v.insert(sys.v.begin(), sys.v_sp.begin(), sys.v_sp.end());
}

void matvec(double* v, sysstruct &sys, double* x)
{
    if (sys.matvecPrecision == 0 && sys.robustMode == 0)
        matvec_sp(v, sys, x);
    else
        matvec_dp(v, sys, x);
}

//void matvec_v3(double* v, sysstruct &sys, double* x)
//{
//    char chn = 'N';
//    double one = 1.0, zero = 0.0;
//
//    Int inc = 1, i, j, nri, pii, pij, rj;
//
//    Int nrows   = sys.numEntities;
//    Int *rowpts = &sys.ent2entStartJ[0];
//    Int *colind = &sys.ent2entJ[0];
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

//void matvec_v4(double* v, sysstruct &sys, double* x)
//{
//    char chn = 'N';
//    double one = 1.0, zero = 0.0;
//
//    Int inc = 1, i, j, nri, pii, pij, rj;
//
//    Int nrows   = sys.numEntities;
//    Int *rowpts = &sys.ent2entStartJ[0];
//    Int *colind = &sys.ent2entJ[0];
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
//void matvec_v5(double* v, sysstruct &sys, double* x)
//{
//    // Note: This implementation showed very slow performance
//    char chn = 'N';
//    double one = 1.0, zero = 0.0;
//
//    Int inc = 1, i, nri, pij = -1, rj, k, kk;
//
//    Int nrows   = sys.numEntities;
//    Int *rowpts = &sys.ent2entStartJ[0];
//    Int *colind = &sys.ent2entJ[0];
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
