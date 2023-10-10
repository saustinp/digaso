#ifndef __MKSHAPG1D
#define __MKSHAPG1D

#include "../master/koornwinder.cpp"

void uniformNodes(double *xi, Int porder)
{
    // Uniform nodes on the interval [0, 1]:
    
    Int numNodes1d = porder + 1;
    
    if (porder == 0)
        xi[0] = 0.5;
    else {
        for (int k = 0; k < numNodes1d; k++)
            xi[k] = ((double) k) / ((double) numNodes1d - 1.0);
    }
}

void extendedChebyshevNodes(double *xi, Int porder)
{
    // Extended Chebyshev nodes on the interval [0, 1]:    
    //printf("Extended Chebyshev nodes not validated yet.\n");
    
    Int numNodes = porder + 1;
    double PI = 3.1415926535897932;
    
    if (porder == 0)
        xi[0] = 0.5;
    else {
        for (int k = 0; k < numNodes; k++) {
            xi[k] = - cos(PI * (2.0 * (double) k + 1.0) / (2.0 * (double) numNodes)) / cos(PI / (2.0 * (double) numNodes));
            xi[k] = 0.5 + 0.5*xi[k];
        }
    }
}

void nodes1d(double *xi, Int porder, Int nodetype)
{
    switch (nodetype) {
        case 0:     // Uniform
            uniformNodes(xi, porder);
            break;
        case 1:     // Extended Chebyshev
            extendedChebyshevNodes(xi, porder);
            break;
        default: {
            cout<<"nodetype not implemented.\n";
            exit(1);
        }
    }
}

void gaussQuad1d(double *x, double *w, Int pgauss)
{
    Int i, ng, lwork, info;
    char chv = 'V', chu = 'U';
    
    if (pgauss % 2)     // pgauss is odd
        ng = (pgauss-1)/2 + 1;
    else                // pgauss is even
        ng = pgauss/2 + 1;
    
    if (ng > 1) {
        lwork = 3*ng-1;
        
        double *beta = new double[ng-1];
        double *T = new double[ng*ng];
        double *L = new double[ng];
        double *work = new double[lwork];

        // 3-term recurrence coefficients for Jacobi matrix:
        for (i = 1; i < ng; i++)
            beta[i-1] = 0.5 / sqrt(1.0 - 1.0/((2.0 * (double) i)*(2.0 * (double) i)));

        // Construct upper triangular part of Jacobi matrix T
        for (i = 0; i < ng*ng; i++)
            T[i] = 0.0;
        for (i = 1; i < ng; i++)
            T[(ng+1)*i-1] = beta[i-1];

        DSYEV(&chv, &chu, &ng, &T[0], &ng, &L[0], &work[0], &lwork, &info);

        // Compute nodes in [0,1] (= Legendre points):
        for (i = 0; i < ng; i++)
            x[i] = 0.5*(L[i] + 1.0);

        // Compute weights in [0 1]:
        for (i = 0; i < ng; i++)
            w[i] = T[i*ng]*T[i*ng];

        delete[] beta; delete[] T;
        delete[] L;
        delete[] work;
    }
    else {
        x[0] = 0.5;
        w[0] = 1.0;
    }
}

void mkshape1d(double* shap, double *plocal, double *pts, Int npoints, Int porder, Int nd, Int numNodes)
{
    // porder: Polynomial order
    // plocal: Node positions. numNodes / nd
    // pts: Points to evaluate shape fucntions and derivatives. npoints / nd
    // shap: shape function and derivatives. numNodes / npoints / nd+1
    
    Int i, j, k, info, inc = 1;
    Int nd1 = nd + 1;
    Int lwork = numNodes;
    char chn = 'N';
    double one = 1.0, zero = 0.0;
    
    Int *ipiv = new Int[numNodes];
    double *work = new double[lwork];
    double *A = new double[numNodes*numNodes];
    double *nf = new double[npoints*numNodes*nd1];

    koornwinder(&nf[0], &pts[0], npoints, porder, nd, 1);           // Orthogonal shape functions
    koornwinder(&A[0], &plocal[0], numNodes, porder, nd, 0);         // Vandermonde matrix
    
    // Divide orthogonal shape functions by the Vandermonde matrix to obtain nodal shape functions:
    DGETRF(&numNodes, &numNodes, &A[0], &numNodes, &ipiv[0], &info);
    double *Ainv = &A[0];
    DGETRI(&numNodes, &Ainv[0], &numNodes, &ipiv[0], &work[0], &lwork, &info);
    for (i = 0; i < nd1; i++) {
        DGEMM(&chn, &chn, &npoints, &numNodes, &numNodes, &one, &nf[i*numNodes*npoints], &npoints, &Ainv[0],
              &numNodes, &zero, &shap[i*numNodes*npoints], &npoints);
        // nf: npoints / numNodes / nd+1
        // Ainv: numNodes / numNodes
        // shap: npoints / numNodes / nd+1
    }
    
    // Permute shap: "npoints / numNodes / nd+1" -> "numNodes / npoints / nd+1"
    double *shap_tmp = new double[numNodes*npoints*(nd+1)];
    for (i = 0; i < nd1; i++)
        for (j = 0; j < numNodes; j++)
            for (k = 0; k < npoints; k++)
                shap_tmp[i*npoints*numNodes+k*numNodes+j] = shap[i*numNodes*npoints+j*npoints+k];
    for (i = 0; i < numNodes*npoints*nd1; i++)
        shap[i] = shap_tmp[i];
    delete[] shap_tmp;
    
    delete[] ipiv; delete[] work;
    delete[] A; delete[] nf;
}

vector<double> mkshapeg1d(Int porder, Int pgauss, Int nodetype)
{
    Int ng;        
    Int nd = 1;
    
    // nodes for 1d shape functions
    vector<double> px(porder+1,0);    
    nodes1d(&px[0], porder, nodetype);    
    
    //  quadrature points in 
    if (pgauss % 2)     // pgauss is odd
        ng = (pgauss-1)/2 + 1;
    else                // pgauss is even
        ng = pgauss/2 + 1;
    vector<double> gx(ng,0);
    vector<double> gw(ng,0);    
    gaussQuad1d(&gx[0], &gw[0], pgauss);
    
    vector<double> shapg((porder+1)*ng*(nd+1),0);    
    mkshape1d(&shapg[0], &px[0], &gx[0], ng, porder, nd, porder+1);
    
    for (int k=0; k<nd+1; k++)
        for (int j=0; j<ng; j++)
            for (int i=0; i<porder+1; i++)
                shapg[k*(porder+1)*ng+j*(porder+1)+i] = shapg[k*(porder+1)*ng+j*(porder+1)+i]*gw[j];    
    
    return shapg;
}

#endif