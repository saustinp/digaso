#ifndef __GMRESLIB
#define __GMRESLIB

// Written by: C. Nguyen & P. Fernandez

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>

using namespace std;

#include "gmres.h"
#include "errormsg.cpp"
#include "pDLAK.cpp"
#include "pSLAK.cpp"
#include "preconditioner.cpp"

static void orthogonalize(double * c, double * q, double * Q, Int m, Int n, Int orthogMethod, 
        double * s, vector<vector<double> >* stmps, Int * convergenceFlag, Int noThreads)
{
    /* AT INPUT:
     * Q (m / n): Set of orthonormal vectors against orthogonalization is to be performed.
     * q (m): Vector to be orthogonalized.
     * m: Number of components of v or q_i.
     * n: Number of vectors in Q.
     * orthogMethod: Orthogonalization method (0: MGS, 1; ICGS, 2: IMGS)
     * s (n): Auxiliary variable for orthogonalization (only required for iterative methods)
     * */

    /* AT OUTPUT:
     * c (n+1): Coefficients of v in the basis of [Q | q].
     * q (m): Orthogonalized version of v.
     * */

    char chn = 'N', cht = 'T';
    double one = 1.0, zero = 0.0, minusone = -1.0;
    double sNorm, qNorm, alpha;
    Int inc = 1, i, j, k, iter, iterPerformed = 0;
    *convergenceFlag = 0;
    
    double relOrthogTol = 1.0e-10;           /* Relative tolerance w.r.t. the norm of q. Only applies for iterative orthogonalization methods (ICGS and IMGS). */
    Int maxIter = 10;                   /* Only applies for iterative orthogonalization methods (ICGS and IMGS). */
    
    if (orthogMethod == 0) {        // MGS method
        for (j=0; j<n; j++) {
            c[j] = DDOT_OpenMP(&m, &Q[j*m], &inc, q, &inc, noThreads);
            alpha = -c[j];
            
            DAXPY_OpenMP(&m, &alpha, &Q[j*m], &inc, q, &inc, noThreads);
        }

        qNorm = DNRM2_OpenMP(&m, q, &inc, noThreads);
        c[n] = qNorm;

        /* Normalize q */
        alpha = 1.0/qNorm;
        DSCAL_OpenMP(&m, &alpha, q, &inc, noThreads);
        
        *convergenceFlag = 1;
    }
    else if (orthogMethod == 1) {       // ICGS method
        // Initialize c to zero:
        for (i=0; i<(n+1); i++)
            c[i] = 0.0;

        /* Iterative orthogonalization */
        iter = 0;
        qNorm = DNRM2_OpenMP(&m, q, &inc, noThreads);
        while (1) {
            /* s = Q^T * q */
            DGEMV_OpenMP(&cht, &m, &n, &one, Q, &m, q, 
                    &inc, &zero, s, &inc, stmps, noThreads);

            /* Check convergence */
            sNorm = 0.0;
            for (i=0; i<n; i++) {
                sNorm += s[i]*s[i];
            }
            sNorm = sqrt(sNorm);            /* Since n is small, better performance may be achieved with for loop instead of BLAS norm routine */

            if (sNorm < relOrthogTol*qNorm) {
                *convergenceFlag = 1;
                break;
            }

            /* Improve q (q = q - Q*s) */
            DGEMV_OpenMP(&chn, &m, &n, &minusone, Q, &m, s,
                    &inc, &one, q, &inc, stmps, noThreads);

            qNorm = DNRM2_OpenMP(&m, q, &inc, noThreads);

            for (i=0; i<n; i++) {
                c[i] += s[i];
            }

            iter += 1;

            if (iter >= maxIter) {
                printf("ICGS orthogonalization did not converge to prescribed tolerance (%g vs. %g) in maxIter = %d.\n", sNorm, relOrthogTol*qNorm, maxIter);
                break;
            }
        }

        iterPerformed = iter;

        c[n] = qNorm;

        /* Normalize q */
        alpha = 1.0/qNorm;
        DSCAL_OpenMP(&m, &alpha, q, &inc, noThreads);
    }
    
    else if (orthogMethod == 2) {       // IMGS method
        // Initialize c to zero:
        for (i=0; i<(n+1); i++)
            c[i] = 0.0;

        /* Iterative orthogonalization */
        iter = 0;
        while (1) {
            iter += 1;

            sNorm = 0.0;
            for (j=0; j<n; j++) {
                s[j] = DDOT_OpenMP(&m, &Q[j*m], &inc, q, &inc, noThreads);

                c[j] += s[j];
                sNorm += s[j]*s[j];

                alpha = -s[j];
                DAXPY_OpenMP(&m, &alpha, &Q[j*m], &inc, q, &inc, noThreads);
            }
            sNorm = sqrt(sNorm);

            qNorm = DNRM2_OpenMP(&m, q, &inc, noThreads);

            /* Check pseudo-convergence (this is a conservative convergence criterion that involves no extra-cost and can avoid checking the expensive true convergence) */
            if (sNorm < relOrthogTol*qNorm) {
                *convergenceFlag = 1;
                break;
            }

            /* Check true convergence (this requires extra work, but it's worth it if convergence is achieved in few iterations due to BLAS2 performance over BLAS1 */
            if (iter <= 10) {
                /* s = Q^T * q */
                DGEMV_OpenMP(&cht, &m, &n, &one, Q, &m, q,
                      &inc, &zero, s, &inc, stmps, noThreads);
                
                sNorm = 0.0;
                for (i=0; i<n; i++) {
                    sNorm += s[i]*s[i];
                }
                sNorm = sqrt(sNorm);

                if (sNorm < relOrthogTol*qNorm) {
                    *convergenceFlag = 1;
                    break;
                }
            }

            if (iter >= maxIter) {
                printf("IMGS orthogonalization did not converge to prescribed tolerance in maxIter = %d\n", maxIter);
                break;
            }
        }

        iterPerformed = iter;

        c[n] = qNorm;

        // Normalize q
        alpha = 1.0/qNorm;
        DSCAL_OpenMP(&m, &alpha, q, &inc, noThreads);
    }
    else {
        printf("orthogMethod flag has invalid value %d.\nOnly 0, 1 and 2 are valid for the serial version of the code.\n", orthogMethod);
        exit(-1);
    }
}

void gmres(sysstruct &sys, Int *convFlag)
{
    // TODO: Implement single/double precision matrix-vector product and orthogonalization in serial version of the code
    if (sys.matvecPrecision == 0)
        error("Single precision matrix-vector product not implemented in serial version of the code.\n");
    if (sys.orthogPrecision == 0)
        error("Single precision orthogonalization not implemented in serial version of the code.\n");
    
    Int iter = 0, cycle = 0, inc = 1, oneInt = 1, info;
    Int i, j, n, lwork, rowsLSP, columnsLSP, flags = -1;
    double alpha, res, beta, one = 1.0, zero = 0.0, nrmb;
    double matvectime = 0, BILU0time = 0, orthogtime = 0;
    char chn = 'N';
    clock_t t, tTotal = clock();
    Int orthConvFlag[1];
    Int * requestCounter;           // Variable only needed for parallel solver
    
    Int preconditionerSide = sys.preconditionerSide;
    
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
    
    Int N = sys.x.size();
    for (i=0; i<N; i++)
        sys.x[i] = 0.0;
    
    /* Save initial residual to Rg_0 */
    for (i=0; i<sys.Rg.size(); i++)
        sys.Rg_0[i] = sys.Rg[i];

    /* Compute preconditioned RHS (only for left preconditioner) */
    if (preconditionerSide == 0)
        applyPreconditioner(sys, &sys.Rg[0], 0, requestCounter);
    nrmb = DNRM2(&N, &sys.Rg[0], &inc);

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
                // Perform matrix-vector multiplication r <- Hg*x
                matvec(&sys.v[0], sys, &sys.x[0]);

                // Apply preconditioner: r <- M*r (=M*Hg*x)
                applyPreconditioner(sys, &sys.v[0], 0, requestCounter);
            }
            else if (preconditionerSide == 1) {
                // Apply preconditioner: x <- M*x
                DCOPY(&N, &sys.x[0], &inc, &sys.Mx[0], &inc);
                applyPreconditioner(sys, &sys.Mx[0], 0, requestCounter);

                // Perform matrix-vector multiplication r <- Hg*x (=Hg*M*x)
                matvec(&sys.v[0], sys, &sys.Mx[0]);
            }

            // r <- M*Rg - r (= M*Rg - M*Hg*x)
            for (i=0; i<N; i++)
                sys.v[i] = sys.Rg[i] - sys.v[i];

            // Compute residual norm
            beta = DNRM2(&N, &sys.v[0], &inc);      /* rNorm should be already available since ||Rg - Hg*M*x_new|| = ||Rg - Hg*M*x_old - Hg*M*e|| = ||r_old - Hg*M*(Q_n*y)|| = || e1*||r_old|| - H^_n * y || */
        }

        // Normalize the residual
        for (i=0; i<N; i++)
            sys.v[i] = sys.v[i]/beta;

        sys.e1[0] = beta;

        for (n=1; n <= sys.restart; n++) {

            if (preconditionerSide == 0) {
                // perform matrix-vector multiplication: r <- Hg*q_n
                t = clock();
                matvec(&sys.v[N*n], sys, &sys.v[N*(n-1)]);
                t = clock() - t;
                matvectime += t;

                // Apply preconditioner: r <- M*r (=M*Hg*q_n)
                t = clock();
                applyPreconditioner(sys, &sys.v[N*n], 0, requestCounter);
                t = clock() - t;
                BILU0time += t;
            }
            else if (preconditionerSide == 1) {
                // Apply preconditioner: Mx <- M*q_n
                DCOPY(&N, &sys.v[N*(n-1)], &inc, &sys.Mv[0], &inc);
                t = clock();
                applyPreconditioner(sys, &sys.Mv[0], 0, requestCounter);
                t = clock() - t;
                BILU0time += t;

                // Perform matrix-vector multiplication r <- Hg*Mx (=Hg*M*q_n)
                t = clock();
                matvec(&sys.v[N*n], sys, &sys.Mv[0]);
                t = clock() - t;
                matvectime += t;
            }
            
            // Arnoldi process (i.e., orthogonalization on the Krylov supspace)
            t = clock();
            orthogonalize(&sys.hh[(n-1)*(sys.restart+1)], &sys.v[N*n], &sys.v[0], N, n, sys.orthogMethod, &sys.s[0], &sys.stmps, &orthConvFlag[0], sys.noThreads);
            if (sys.hh[(n-1)*(sys.restart+1)+n] == 0.0) {
                printf("\n\nWARNING: A breakdown in the GMRES method occurred.\n\n");
                flags = 2;
                break;
            }
            else if (orthConvFlag[0] == 0) {
                flags = 2;
                break;
            }
            t = clock() - t;
            orthogtime += t;

            rowsLSP = n+1;
            columnsLSP = n;
            for (j=0; j<columnsLSP; j++)
                for (i=0; i<rowsLSP; i++)
                    sys.hhls[j*rowsLSP+i] = sys.hh[j*(sys.restart+1)+i];

            for (i=0; i<rowsLSP; i++)
                sys.y[i] = sys.e1[i];

            // solve the least-squares problem (reduced system)
            //y = hh(1:j+1,1:j)\(beta*e1(1:j+1));
            lwork = rowsLSP + columnsLSP;
            DGELS(&chn, &rowsLSP, &columnsLSP, &oneInt, &sys.hhls[0], &rowsLSP, &sys.y[0], &rowsLSP, &sys.workls[0], &lwork, &info);      /* Note: At output, hhls is overwritten by details of its QR factorization. */

            for (j=0; j<columnsLSP; j++)
                for (i=0; i<rowsLSP; i++)
                    sys.hhls[j*rowsLSP+i] = sys.hh[j*(sys.restart+1)+i];

            DGEMV(&chn, &rowsLSP, &columnsLSP, &one, &sys.hhls[0], &rowsLSP, &sys.y[0],
                  &inc, &zero, &sys.hy[0], &inc);

            // compute the preconditioned residual norm
            for (i=0; i<rowsLSP; i++)
                sys.hy[i] = sys.e1[i] - sys.hy[i];
            res = DNRM2(&rowsLSP, &sys.hy[0], &inc);

            sys.rev[iter] = res;
            iter = iter + 1;

            if (sys.print >= 3) {
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

        if (flags == 0 || flags == 1) {
        }
        else {
            n -= 1;
        }

        // Compute the solution:
        DGEMV_OpenMP(&chn, &N, &n, &one, &sys.v[0], &N, &sys.y[0],
                  &inc, &one, &sys.x[0], &inc, &sys.stmps, sys.noThreads);

        cycle = cycle + 1;

        // Stop if converging, reaching the maximum iteration or failure (i.e. breakdown, orthogonalization does not converge)
        if (flags == 0 || flags == 1 || flags == 2) {
            break;
        }
    }

    /* Apply right preconditioner to recover final solution */
    if (preconditionerSide == 1)
        applyPreconditioner(sys, &sys.x[0], 0, requestCounter);

    // Report GMRES outcome:
    if (flags == 0) {
        *convFlag = 1;
        printf("\nGMRES(%d) converges at iteration %d to a solution with relative residual %g\n", sys.restart, iter, res/nrmb);
    }
    else if (flags == 1) {
        *convFlag = 0;
        printf("\nGMRES(%d) did not converge in maximum number of iterations (%d). Relative residual %g. Prescribed tolerance %g\n", sys.restart, iter, res/nrmb, GMREStol);
    }
    else if (flags == 2)
        *convFlag = -1;

    Int bytesPerNumber_matvec, bytesPerNumber_precSolve;
    if (sys.matvecPrecision == 0 && sys.robustMode == 0)
        bytesPerNumber_matvec = 4;
    else
        bytesPerNumber_matvec = 8;
    if (sys.precPrecision == 0 && sys.robustMode == 0)
        bytesPerNumber_precSolve = 4;
    else
        bytesPerNumber_precSolve = 8;
    if (sys.print >= 2) {
        printf("Matrix-vector time: %g ms || Preconditioner solve time: %g ms || Orthogonalization time: %g ms\n", (matvectime/CLOCKS_PER_SEC)*1.0e3, (BILU0time/CLOCKS_PER_SEC)*1.0e3, (orthogtime/CLOCKS_PER_SEC)*1.0e3);
        printf("Preconditioner solve performance: %g GBytes/s \n", (( (double) (bytesPerNumber_precSolve * iter * (long) sys.blkSize * (long) sys.blkSize * (long) sys.ent2entStart[sys.numEntities]))/1.0e9)/(BILU0time/CLOCKS_PER_SEC));
        printf("Matrix-vector product performance: %g GBytes/s \n", (( (double) (bytesPerNumber_matvec * iter * (long) sys.blkSize * (long) sys.blkSize*sys.ent2entStart[sys.numEntities]))/1.0e9)/(matvectime/CLOCKS_PER_SEC));
    }
}

#endif
