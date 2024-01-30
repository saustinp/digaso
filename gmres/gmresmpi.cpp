#ifndef __GMRESMPI
#define __GMRESMPI

// Written by: C. Nguyen & P. Fernandez
#include "sendrecvMPI.cpp"
#include "pDLAK.cpp"
#include "pSLAK.cpp"
#include "preconditioner.cpp"

void applyMPIpreconditioner(sysstruct &sys, double * r, Int preconditioner, Int sendInputVector, Int sendOutputVector, Int * requestCounter, double* precondtimes)
{
    clock_t t;
    
    // Send input vector (only required for RAS preconditioner)
    t = clock();
    if (preconditioner == 0 && sendInputVector == 1)
        sendvector(sys, &r[0], requestCounter);
    precondtimes[0] += clock() - t;

    // Apply preconditioner
    t = clock();
    applyPreconditioner(sys, &r[0], preconditioner, requestCounter);
    precondtimes[1] += clock() - t;

    // Send output vector (if necessary). Reception will be done is other part of the code
    if (sendOutputVector == 1)
        sendvector(sys, &r[0], requestCounter);
}

void orthogonalizempi_dp(double* c, double* q, double* Q, Int m, Int n, Int orthogMethod, 
        double* s, double* stmp, vector<vector<double> >* stmps, Int* convergenceFlag, Int* iterPerformedPt,
        Int my_rank, Int noThreads)
{
    char chn = 'N', cht = 'T';
    double one = 1.0, zero = 0.0, minusone = -1.0;
    double sNorm, qNorm, alpha;
    Int inc = 1, i, j, k, iter;
    *convergenceFlag = 0;
    
    double relOrthogTol = 1.0e-6;           /* Relative tolerance w.r.t. the norm of q. Only applies for iterative orthogonalization methods. */
    Int maxIter = 10;                       /* Only applies for iterative orthogonalization methods. */
    
    // Orthogonalize:
    if (orthogMethod == 0) {        /* MGS method */
        for (j=0; j<n; j++) {
            c[j] = DDOTMPI(m, &Q[j*m], inc, q, inc, noThreads);
            alpha = -c[j];
            DAXPY_OpenMP(&m, &alpha, &Q[j*m], &inc, q, &inc, noThreads);
        }

        qNorm = sqrt(DDOTMPI(m, q, inc, q, inc, noThreads));        

        c[n] = qNorm;

        // Normalize q
        alpha = 1.0 / qNorm;
        DSCAL_OpenMP(&m, &alpha, q, &inc, noThreads);
        
        *convergenceFlag = 1;
    }
    else if (orthogMethod == 1) {       /* ICGS method */
        // Initialize c to zero:
        for (i=0; i<(n+1); i++)
            c[i] = 0.0;
        
        // Iterative orthogonalization
        iter = 0;
        while (1) {
            /* s = Q^T * q */
            DGEMVMPI(cht, m, n+1, one, Q, m, q, inc, zero, s, inc, &stmp[0], stmps, noThreads);
            
            // Check convergence:
            sNorm = 0.0;
            for (i=0; i<n; i++) {
                sNorm += s[i]*s[i];
            }
            sNorm = sqrt(sNorm);            /* Since n is small, better performance may be achieved with for loop instead of BLAS norm routine */
            qNorm = sqrt(s[n]);

            if (sNorm < relOrthogTol*qNorm) {
                *convergenceFlag = 1;
                break;
            }
            
            if (iter >= maxIter) {
                if (my_rank == 0) {
                    printf("ICGS orthogonalization did not converge to prescribed tolerance (%g vs. %g) in maxIter = %d.\n", sNorm, relOrthogTol*qNorm, maxIter);
                }
                break;
            }

            // Improve q (q = q - Q*s)
            DGEMV_OpenMP(&chn, &m, &n, &minusone, Q, &m, s,
                  &inc, &one, q, &inc, stmps, noThreads);
            
            for (i=0; i<n; i++) {
                c[i] += s[i];
            }
            iter += 1;
        }
        * iterPerformedPt = iter;
        
        c[n] = qNorm;

        /* Normalize q */
        alpha = 1.0 / qNorm;
        DSCAL_OpenMP(&m, &alpha, q, &inc, noThreads);
    }
    else if (orthogMethod == 2) {
        if (my_rank == 0)
            printf("IMGS not implemented in the MPI code yet.\n", orthogMethod);
        error("\n");
    }
    else if (orthogMethod == 3) {       /* ICGS method with no convergence guarantee */
        maxIter = *iterPerformedPt;
        
        // Initialize c to zero:
        for (i=0; i<(n+1); i++)
            c[i] = 0.0;

        // Iterative orthogonalization
        iter = 0;
        while (1) {
            /* s = Q^T * q */
            DGEMVMPI(cht, m, n+1, one, Q, m, q, inc, zero, s, inc, &stmp[0], stmps, noThreads);
            
            // Check convergence
            sNorm = 0.0;
            for (i=0; i<n; i++) {
                sNorm += s[i]*s[i];
            }
            sNorm = sqrt(sNorm);            /* Since n is small, better performance may be achieved with for loop instead of BLAS norm routine */
            qNorm = sqrt(s[n]);

            if (sNorm < relOrthogTol*qNorm) {
                *convergenceFlag = 1;
                break;
            }
            
            // Improve q (q = q - Q*s)
            DGEMV_OpenMP(&chn, &m, &n, &minusone, Q, &m, s,
                  &inc, &one, q, &inc, stmps, noThreads);
            
            for (i=0; i<n; i++) {
                c[i] += s[i];
            }
            iter += 1;

            // Check maximum number of iterations
            if (iter >= maxIter) {
                qNorm = sqrt(DDOTMPI(m, q, inc, q, inc, noThreads));
                *convergenceFlag = 1;
                break;
            }
        }

        if (iter > 1)
            * iterPerformedPt = iter;
        else
            * iterPerformedPt = 1;

        c[n] = qNorm;

        // Normalize q
        alpha = 1.0 / qNorm;
        DSCAL_OpenMP(&m, &alpha, q, &inc, noThreads);
    }
    else {
        if (my_rank == 0)
            printf("orthogMethod flag has invalid value %d.\nOnly 0, 1 and 3 are valid for the MPI version of the code.\n", orthogMethod);
        error("\n");
    }
}

void orthogonalizempi_sp(double* c, double* q, double* Q, float* q_sp, float* Q_sp, Int m, Int n, Int orthogMethod, 
        float* s_sp, float* stmp_sp, vector<vector<float> >* stmps_sp, Int* convergenceFlag, Int* iterPerformedPt, Int my_rank, Int noThreads)
{
    char chn = 'N', cht = 'T';
    float one = 1.0, zero = 0.0, minusone = -1.0;
    float sNorm, qNorm, alpha;
    Int inc = 1, i, j, k, iter;
    *convergenceFlag = 0;
    
    float relOrthogTol = 1.0e-6;           /* Relative tolerance w.r.t. the norm of q. Only applies for iterative orthogonalization methods. */
    Int maxIter = 10;                       /* Only applies for iterative orthogonalization methods. */
    
    // Convert Q, q to single precision:
    for (i = (n-1)*m; i < n*m; i++)
        Q_sp[i] = (float) Q[i];
    for (i = 0; i < m; i++)
        q_sp[i] = (float) q[i];
    
    // Orthogonalize:
    if (orthogMethod == 0) {        /* MGS method */    
        for (j=0; j<n; j++) {
            alpha = - SDOTMPI(m, &Q_sp[j*m], inc, q_sp, inc, noThreads);
            c[j] = (double) - alpha;
            SAXPY_OpenMP(&m, &alpha, &Q_sp[j*m], &inc, q_sp, &inc, noThreads);
        }

        qNorm = sqrt(SDOTMPI(m, q_sp, inc, q_sp, inc, noThreads));        

        c[n] = (double) qNorm;

        // Normalize q
        alpha = 1.0 / qNorm;
        SSCAL_OpenMP(&m, &alpha, q_sp, &inc, noThreads);
        
        *convergenceFlag = 1;
    }
    else if (orthogMethod == 1) {       /* ICGS method */
        // Initialize c to zero:
        for (i=0; i<(n+1); i++)
            c[i] = 0.0;
        
        // Iterative orthogonalization
        iter = 0;
        while (1) {
            /* s = Q^T * q */
            SGEMVMPI(cht, m, n+1, one, Q_sp, m, q_sp, inc, zero, s_sp, inc, stmp_sp, stmps_sp, noThreads);
            
            // Check convergence:
            sNorm = 0.0;
            for (i=0; i<n; i++) {
                sNorm += s_sp[i]*s_sp[i];
            }
            sNorm = sqrt(sNorm);            /* Since n is small, better performance may be achieved with for loop instead of BLAS norm routine */
            qNorm = sqrt(s_sp[n]);

            if (sNorm < relOrthogTol*qNorm) {
                *convergenceFlag = 1;
                break;
            }
            
            if (iter >= maxIter) {
                if (my_rank == 0) {
                    printf("ICGS orthogonalization did not converge to prescribed tolerance (%g vs. %g) in maxIter = %d.\n", sNorm, relOrthogTol*qNorm, maxIter);
                }
                break;
            }

            // Improve q (q = q - Q*s)
            SGEMV_OpenMP(&chn, &m, &n, &minusone, Q_sp, &m, s_sp,
                  &inc, &one, q_sp, &inc, stmps_sp, noThreads);
            
            for (i=0; i<n; i++) {
                c[i] += (double) s_sp[i];
            }
            iter += 1;
        }

        * iterPerformedPt = iter;
        
        c[n] = (double) qNorm;

        /* Normalize q */
        alpha = 1.0 / qNorm;
        SSCAL_OpenMP(&m, &alpha, q_sp, &inc, noThreads);
    }
    else if (orthogMethod == 2) {
        if (my_rank == 0)
            printf("IMGS not implemented in the MPI code yet.\n", orthogMethod);
        error("\n");
    }
    else if (orthogMethod == 3) {       /* ICGS method with no convergence guarantee */
        maxIter = * iterPerformedPt;

        // Initialize c to zero:
        for (i=0; i<(n+1); i++)
            c[i] = 0.0;

        // Iterative orthogonalization
        iter = 0;
        while (1) {
            /* s = Q^T * q */
            SGEMVMPI(cht, m, n+1, one, Q_sp, m, q_sp, inc, zero, s_sp, inc, stmp_sp, stmps_sp, noThreads);
            
            // Check convergence
            sNorm = 0.0;
            for (i=0; i<n; i++) {
                sNorm += s_sp[i]*s_sp[i];
            }
            sNorm = sqrt(sNorm);            /* Since n is small, better performance may be achieved with for loop instead of BLAS norm routine */
            qNorm = sqrt(s_sp[n]);

            if (sNorm < relOrthogTol*qNorm) {
                *convergenceFlag = 1;
                break;
            }

            // Improve q (q = q - Q*s)
            SGEMV_OpenMP(&chn, &m, &n, &minusone, Q_sp, &m, s_sp,
                  &inc, &one, q_sp, &inc, stmps_sp, noThreads);
            
            for (i=0; i<n; i++) {
                c[i] += (double) s_sp[i];
            }
            iter += 1;

            // Check maximum number of iterations
            if (iter >= maxIter) {
                qNorm = sqrt(SDOTMPI(m, q_sp, inc, q_sp, inc, noThreads));
                *convergenceFlag = 1;
                break;
            }
        }

        if (iter > 1)
            * iterPerformedPt = iter;
        else
            * iterPerformedPt = 1;

        c[n] = (double) qNorm;

        // Normalize q:
        alpha = 1.0 / qNorm;
        SSCAL_OpenMP(&m, &alpha, q_sp, &inc, noThreads);
    }
    else {
        if (my_rank == 0)
            printf("orthogMethod flag has invalid value %d.\nOnly 0, 1 and 3 are valid for the MPI version of the code.\n", orthogMethod);
        error("\n");
    }
    
    // Convert q to double precision:
    for (i = 0; i < m; i++)
        q[i] = (double) q_sp[i];
}


void orthogonalizempi(double* c, double* q, double* Q, float* q_sp, float* Q_sp, 
        Int m, Int n, Int orthogMethod, double* s, double* stmp, vector<vector<double> >* stmps,
        float* s_sp, float* stmp_sp, vector<vector<float> >* stmps_sp, Int* convergenceFlag, Int* iterPerformedPt, sysstruct &sys)
{
    /* AT INPUT:
     * Q (m / n): Set of orthonormal vectors against orthogonalization is to be performed.
     * q (m): Vector to be orthogonalized.
     * m: Number of components of v or q_i.
     * n: Number of vectors in Q.
     * orthogMethod: Orthogonalization method (0: MGS, 1; ICGS, 2: IMGS, 2: ICGS with no convergece guarantee)
     * s (n): Auxiliary variable for orthogonalization (only required for iterative methods)
     * */

    /* AT OUTPUT:
     * c (n+1): Coefficients of v in the basis of [Q | q].
     * q (m): Orthogonalized version of v.
     * */

    if (sys.orthogPrecision == 0 && sys.robustMode == 0)
        orthogonalizempi_sp(c, q, Q, q_sp, Q_sp, m, n, orthogMethod, s_sp, stmp_sp, stmps_sp,
                convergenceFlag, iterPerformedPt, sys.my_rank, sys.noThreads);
    else
        orthogonalizempi_dp(c, q, Q, m, n, orthogMethod, s, stmp, stmps,
                convergenceFlag, iterPerformedPt, sys.my_rank, sys.noThreads);
}

void gmresmpi(sysstruct &sys, Int *convFlag)
{
    Int iter = 0, cycle = 0, inc = 1, oneInt = 1, info;
    Int i, j, n, lwork, rowsLSP, columnsLSP, flags = -1;
    Int orthIterPerformed;
    Int numGMRESiterToCheckICGSconv = 10;       // Number of GMRES iterations after which it is actually verified that ICGS converges. Only applies for orthogMethod == 3.
    double alpha, res, beta, one = 1.0, zero = 0.0, nrmb;
    double matvectime = 0, BILU0time = 0, orthogtime = 0;
    clock_t t, tTotal = clock();
    char chn = 'N';
    Int orthConvFlag[1];

    Int N = sys.BJ_nrows*sys.blkSize;
    Int preconditionerSide = sys.preconditionerSide;
    Int preconditioner = sys.preconditioner;
    Int sendOutputVector, sendInputVector;
    Int requestCounter[1];
    
    //error("before gmres iter -1");
    double matvectimes[3];
    for (i = 0; i < 3; i++) {
        matvectimes[i] = 0.0;
    }
    double precondtimes[2];
    for (i = 0; i < 2; i++) {
        precondtimes[i] = 0.0;
    }
    
    // Decide GMRES tolerance:
    double GMREStol;
    double securityFactor = 50.0;       // Security factor to ensure that (if possible) non-linear convergence is achieved in the current Newton iteration.
    if (sys.linearProblem == 0 && sys.robustMode == 1) {
        GMREStol = min(sys.tol,1.0e-9);
        sys.adaptGMREStol = 1;
    }
    else if (sys.linearProblem == 0 && sys.robustMode == 0 && sys.adaptiveGMREStol == 1 && sys.adaptGMREStol == 1 && ((sys.NewtonTol / (securityFactor * sys.rNorm)) > sys.tol)) {
        GMREStol = sys.NewtonTol / (securityFactor * sys.rNorm);
        sys.adaptGMREStol = 1;
    }
    else if (sys.linearProblem == 0 && sys.adaptiveGMREStol == 1 && sys.adaptGMREStol == 0) {
        // Adaptive tolerance did not allow converging Newton in previous iteration. Reduce linear tolerance to ensure descent direction.
        GMREStol = min(sys.tol,1.0e-9);
        sys.adaptGMREStol = 1;
    }
    else {
        GMREStol = sys.tol;
        sys.adaptGMREStol = 0;
    }

    
    //error("before gmres iter 0");
    
    /* Initial guess for GMRES */
    for (i=0; i<sys.x.size(); i++)
        sys.x[i] = 0.0;
    
    /* Save initial residual to Rg_0 */
    for (i=0; i<sys.Rg.size(); i++)
        sys.Rg_0[i] = sys.Rg[i];
        
    
    /* Compute preconditioned RHS (only for left preconditioner) */
    if (preconditionerSide == 0) {
        t = clock();
        sendInputVector = 1;
        sendOutputVector = 0;
        applyMPIpreconditioner(sys, &sys.Rg[0], preconditioner, sendInputVector, sendOutputVector, requestCounter, &precondtimes[0]);     // TODO: Need to allocate more memory for Rg (??)
        BILU0time += clock() - t;
    }
    
    //error("before gmres iter 1");
    
    nrmb = sqrt(DDOTMPI(N, &sys.Rg[0], inc, &sys.Rg[0], inc, sys.noThreads));
    
    //error("before gmres iter");
    while (iter < sys.maxiter) {     
        
        if (iter == 0) {
            // Residual is already available
            DCOPY(&N, &sys.Rg[0], &inc, &sys.v[0], &inc);
            beta = nrmb;
            if (beta == 0.0) {
                // RHS is zero and initialization is solution
                flags = 0;
                res = 0.0;
                break;
            }
        }
        else {
            if (preconditionerSide == 0) {
                // Perform matrix-vector multiplication Mx <- (Hg+Kg)*x
                sendInputVector = 1;       // sys.x was just computed (no intermediate operations before reaching this line)
                if (preconditioner == 0)
                    sendOutputVector = 1;
                else if (preconditioner == -1 || preconditioner == 1 || preconditioner == 2)
                    sendOutputVector = 0;
                matvecMPI(&sys.Mx[0], sys, &sys.x[0], sendInputVector, sendOutputVector, requestCounter, &matvectimes[0]);
                
                // Apply preconditioner: Mx <- M*r (=M*(Hg+Kg)*x)
                t = clock();
                sendInputVector = 0;
                sendOutputVector = 0;
                applyMPIpreconditioner(sys, &sys.Mx[0], preconditioner, sendInputVector, sendOutputVector, requestCounter, &precondtimes[0]);
                BILU0time += clock() - t;
                
                DCOPY(&N, &sys.Mx[0], &inc, &sys.v[0], &inc);
            }
            else if (preconditionerSide == 1) {                
                DCOPY(&N, &sys.x[0], &inc, &sys.Mx[0], &inc);
                
                // Apply preconditioner: Mx <- M*x
                t = clock();
                sendInputVector = 1;       // sys.x (and hence sys.Mx) was just computed (no intermediate operations before reaching this line)
                sendOutputVector = 1;
                applyMPIpreconditioner(sys, &sys.Mx[0], preconditioner, sendInputVector, sendOutputVector, requestCounter, &precondtimes[0]);
                BILU0time += clock() - t;
                
                // Perform matrix-vector multiplication r <- (Hg+Kg)*Mx (=(Hg+Kg)*M*x)
                sendInputVector = 0;
                sendOutputVector = 0;
                matvecMPI(&sys.v[0], sys, &sys.Mx[0], sendInputVector, sendOutputVector, requestCounter, &matvectimes[0]);
            }
            
            // r <- M*Rg - r (= M*Rg - M*Hg*x)
            for (i=0; i<N; i++)
                sys.v[i] = sys.Rg[i] - sys.v[i];

            // Compute residual norm
            beta = sqrt(DDOTMPI(N, &sys.v[0], inc, &sys.v[0], inc, sys.noThreads));
        }

        // Normalize the residual
        for (i=0; i<N; i++)
            sys.v[i] = sys.v[i]/beta;

        // Send q1 = &sys.v[0] among processors
        sendvector(sys, &sys.v[0], requestCounter);

        sys.e1[0] = beta;

        
        for (n=1; n <= sys.restart; n++) {

            if (preconditionerSide == 0) {  
                // Perform matrix-vector multiplication r <- (Hg+Kg)*x
//                 if (sys.my_rank == 1) {
//                  printf("1st GMRES %d with relative residual %g\n", iter, res/nrmb);
//                  cout<<sys.v.size()<<endl;
//                  cout<<N*(n-1)<<endl;
//                 }
                t = clock();
                sendInputVector = 0;
                if (preconditioner == 0)
                    sendOutputVector = 1;
                else if (preconditioner == -1 || preconditioner == 1 || preconditioner == 2)
                    sendOutputVector = 0;
                
                matvecMPI(&sys.Mx[0], sys, &sys.v[N*(n-1)], sendInputVector, sendOutputVector, requestCounter, &matvectimes[0]);
                matvectime += clock() - t;
                
                MPI_Barrier(MPI_COMM_WORLD);          
                
//                 if (sys.my_rank == 1)
//                  printf("2nd GMRES %d with relative residual %g\n", iter, res/nrmb);
                
                // Apply preconditioner: r <- M*r (=M*(Hg+Kg)*x)
                t = clock();
                sendInputVector = 0;
                sendOutputVector = 0;
                applyMPIpreconditioner(sys, &sys.Mx[0], preconditioner, sendInputVector, sendOutputVector, requestCounter, &precondtimes[0]);
                BILU0time += clock() - t;
                
//                 if (sys.my_rank == 1)
//                  printf("3rd GMRES %d with relative residual %g\n", iter, res/nrmb);
//                 
                DCOPY(&N, &sys.Mx[0], &inc, &sys.v[N*n], &inc);     // TODO: Thread this line and other similar lines in GMRES
                
//                 if (sys.my_rank == 1)
//                  printf("4th GMRES %d with relative residual %g\n", iter, res/nrmb);                                                
                                
                MPI_Barrier(MPI_COMM_WORLD);    
            }
            else if (preconditionerSide == 1) {                
                DCOPY(&N, &sys.v[N*(n-1)], &inc, &sys.Mx[0], &inc);

                // Apply preconditioner: x <- M*x
                t = clock();
                sendInputVector = 0;
                sendOutputVector = 1;
                applyMPIpreconditioner(sys, &sys.Mx[0], preconditioner, sendInputVector, sendOutputVector, requestCounter, &precondtimes[0]);
                BILU0time += clock() - t;
                
                // Perform matrix-vector multiplication r <- Hg*x (=Hg*M*x)
                t = clock();
                sendInputVector = 0;
                sendOutputVector = 0;
                matvecMPI(&sys.v[N*n], sys, &sys.Mx[0], sendInputVector, sendOutputVector, requestCounter, &matvectimes[0]);
                matvectime += clock() - t;
            }            
    
//             if (sys.my_rank == 1)
//              printf("5th GMRES %d with relative residual %g\n", iter, res/nrmb);                                                
            
            // Arnoldi process (orthogonalization in the Krylov subspace)
            t = clock();
            if (sys.orthogMethod == 3 && ((n-1) % numGMRESiterToCheckICGSconv) == 0) {       // Set a large value of maximum ICGS iterations
                orthIterPerformed = 10;
            }
            orthogonalizempi(&sys.hh[(n-1)*(sys.restart+1)], &sys.v[N*n], &sys.v[0], &sys.v_sp[N*n], &sys.v_sp[0], N, n, sys.orthogMethod, 
                    &sys.s[0], &sys.stmp[0], &sys.stmps, &sys.s_sp[0], &sys.stmp_sp[0], &sys.stmps_sp, &orthConvFlag[0], &orthIterPerformed, sys);
            
//             if (sys.my_rank == 1)
//              printf("6th GMRES %d with relative residual %g\n", iter, res/nrmb);                                                
            
            if (sys.hh[(n-1)*(sys.restart+1)+n] == 0.0) {
                if (sys.my_rank == 1)
                    printf("\n\nWARNING: A breakdown in the GMRES method occurred.\n\n");
                flags = 2;
                break;
            }
            else if (orthConvFlag[0] == 0) {
                if (sys.my_rank == 1)
                    printf("\n\nWARNING: Orthogonalization in GMRES did not converge.\n\n");
                flags = 2;
                break;
            }
            t = clock() - t;
            orthogtime += t;
            
            // Send q_n = &sys.v[N*n] among processors
            sendvector(sys, &sys.v[N*n], requestCounter);

            rowsLSP = n+1;
            columnsLSP = n;
            for (j=0; j<columnsLSP; j++)
                for (i=0; i<rowsLSP; i++) 
                    sys.hhls[j*rowsLSP+i] = sys.hh[j*(sys.restart+1)+i];                
            
            for (i=0; i<rowsLSP; i++)
                sys.y[i] = sys.e1[i];
            
            // solve the least-squares problem (reduced system)
            lwork = rowsLSP + columnsLSP;
            DGELS(&chn, &rowsLSP, &columnsLSP, &oneInt, &sys.hhls[0], &rowsLSP, &sys.y[0], &rowsLSP, &sys.workls[0], &lwork, &info);      /* Note: At output, hhls is overwritten by details of its QR factorization. */
            
//             if (sys.my_rank == 1)
//              printf("7th GMRES %d with relative residual %g\n", iter, res/nrmb);                                                
            
            for (j=0; j<columnsLSP; j++)
                for (i=0; i<rowsLSP; i++) 
                    sys.hhls[j*rowsLSP+i] = sys.hh[j*(sys.restart+1)+i];                
            
            DGEMV(&chn, &rowsLSP, &columnsLSP, &one, &sys.hhls[0], &rowsLSP, &sys.y[0],
                  &inc, &zero, &sys.hy[0], &inc);
            
            // Compute the preconditioned residual norm
            for (i=0; i<rowsLSP; i++) 
                sys.hy[i] = sys.e1[i] - sys.hy[i];
            
            res = DNRM2(&rowsLSP, &sys.hy[0], &inc);
            
            sys.rev[iter] = res;
            iter = iter + 1;
            
            if (sys.print >= 3) {
                if (sys.my_rank == 0)
                    printf("GMRES Iter No. %d with relative residual %g\n", iter, res/nrmb);
            }

            // Set flag=0 if convergence
            if (res/nrmb < GMREStol) {
                flags = 0;
                break;
            }
            
            // Set flag=1 if reaching maximum iteration
            if (iter >= sys.maxiter) {
                flags = 1;
                break;
            }
        }
        
        //error("here");

        if (flags == 0 || flags == 1) {
        }
        else {
            n -= 1;
        }
        
        // Compute the solution
        DGEMV_OpenMP(&chn, &N, &n, &one, &sys.v[0], &N, &sys.y[0],
                  &inc, &one, &sys.x[0], &inc, &sys.stmps, sys.noThreads);
        
        cycle = cycle + 1;
        
        // Receive q_n = &sys.v[N*n]. TODO: The following receive doesn't seem to be required for the code to work, but it's kept just in case
        recvvector(sys, &sys.v[N*n], requestCounter);
        
        // Stop if converging, reaching the maximum iteration or failure (i.e. breakdown, orthogonalization does not converge)
        if (flags == 0 || flags == 1 || flags == 2) {
            break;
        }
    }
    
    /* Apply right preconditioner to recover final solution */
    if (preconditionerSide == 1) {
        t = clock();
        sendInputVector = 1;       // sys.x was just computed (no intermediate operations before reaching this line)
        sendOutputVector = 0;
        applyMPIpreconditioner(sys, &sys.x[0], preconditioner, sendInputVector, sendOutputVector, requestCounter, &precondtimes[0]);
        BILU0time += clock() - t;
    }
    
    // Report GMRES outcome:
    if (flags == 0) {
        *convFlag = 1;
        if (sys.my_rank == 0)
            printf("\nGMRES(%d) converges at iteration %d to a solution with relative residual %g\n", sys.restart, iter, res/nrmb);
    }
    else if (flags == 1) {
        *convFlag = 0;
        if (sys.my_rank == 0)
            printf("\nGMRES(%d) did not converge in maximum number of iterations (%d). Relative residual %g. Prescribed tolerance %g\n", sys.restart, iter, res/nrmb, GMREStol);
    }
    else if (flags == 2) {
        *convFlag = -1;
        if (sys.my_rank == 0)
            printf("\nGMRES failure.\n");
    }

    Int bytesPerNumber_matvec, bytesPerNumber_precSolve;
    if (sys.matvecPrecision == 0 && sys.robustMode == 0)
        bytesPerNumber_matvec = 4;
    else
        bytesPerNumber_matvec = 8;
    if (sys.precPrecision == 0 && sys.robustMode == 0)
        bytesPerNumber_precSolve = 4;
    else
        bytesPerNumber_precSolve = 8;
    if (sys.print >= 2 && sys.my_rank == 0) {
        Int blocksInPreconditioner;

        printf("\nBREAKDOWN OF GMRES CPU TIME:\n");
        printf("\n1. Matrix-vector time: %g ms\n", (matvectime/CLOCKS_PER_SEC)*1.0e3);
        for (i = 0; i < 3; i++) {
            printf("1.%d: %g ms\n", i, (matvectimes[i]/CLOCKS_PER_SEC)*1.0e3);
        }
        printf("Matrix-vector product performance: %g GBytes/s \n", (bytesPerNumber_matvec * iter * (long) sys.blkSize * (long) sys.blkSize*sys.ent2entStart[sys.BJ_nrows]/1.0e9)/(matvectime/CLOCKS_PER_SEC));

        printf("\n2. Preconditioner solve time: %g ms\n", (BILU0time/CLOCKS_PER_SEC)*1.0e3);
        for (i = 0; i < 2; i++) {
            printf("2.%d: %g ms\n", i, (precondtimes[i]/CLOCKS_PER_SEC)*1.0e3);
        }
        
        if (preconditioner == -1)
            blocksInPreconditioner = 0;
        else if (preconditioner == 0)
            blocksInPreconditioner = sys.nblksP;
        else if (preconditioner == 1)
            blocksInPreconditioner = sys.ent2entStart[sys.BJ_nrows];
        else if (preconditioner == 2)
            blocksInPreconditioner = sys.BJ_nrows;
        if (preconditioner != -1)
            printf("Preconditioner performance: %g GBytes/s \n", (bytesPerNumber_precSolve * iter * (long) sys.blkSize * (long) sys.blkSize * (long) blocksInPreconditioner/1.0e9)/(BILU0time/(CLOCKS_PER_SEC)));
        
        printf("\n3. Orthogonalization time: %g ms\n", orthogtime*1.0e3/CLOCKS_PER_SEC);
    }
}

#endif
