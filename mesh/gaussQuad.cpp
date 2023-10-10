#ifndef __GAUSSQUAD
#define __GAUSSQUAD

#include "gaussQuadTri.cpp"
#include "gaussQuadTet.cpp"

// Written by: C. Nguyen & P. Fernandez

void gaussQuad1D(vector<double> &x, vector<double> &w, Int pgauss)
{
    Int i, ng, lwork, info;
    char chv = 'V', chu = 'U';
    
    if (pgauss % 2)     // pgauss is odd
        ng = (pgauss-1)/2 + 1;
    else                // pgauss is even
        ng = pgauss/2 + 1;
    
    x.resize(ng,0);
    w.resize(ng,0);    
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

void gaussQuadTensor(vector<double> &xg, vector<double> &wg, vector<Int> pgauss)
{
    // x: ng / nd
    // w: ng
    
    Int i, j, k, ng, ng2d, dim = pgauss.size();
    vector<Int> ng1d(dim,0);
    vector< vector<double> > x1d(dim, vector<double>());
    vector< vector<double> > w1d(dim, vector<double>());
            
    for (i=0; i<dim; i++) {
        gaussQuad1D(x1d[i], w1d[i], pgauss[i]);
        ng1d[i] = w1d.size();
    }
    
    switch (dim) {
        case 1:
            xg = x1d[0];
            wg = w1d[0];
            break;
        case 2:
            ng = ng1d[0]*ng1d[1];
            xg.resize(ng*2);
            wg.resize(ng);
            for (i = 0; i < ng1d[1]; i++) 
                for (j = 0; j < ng1d[0]; j++) {
                    xg[0*ng+i*ng1d[0]+j] = x1d[0][j];
                    xg[1*ng+i*ng1d[0]+j] = x1d[1][i];
                    wg[i*ng1d[0]+j] = w1d[1][i]*w1d[0][j];
                }                    
            break;
        case 3:
            ng = ng1d[0]*ng1d[1]*ng1d[2];
            xg.resize(ng*3);
            wg.resize(ng);
            ng2d = ng1d[0]*ng1d[1];
            for (i = 0; i < ng1d[2]; i++) 
                for (j = 0; j < ng1d[1]; j++) 
                    for (k = 0; k < ng1d[0]; k++) {
                        xg[0*ng+i*ng2d+j*ng1d[0]+k] = x1d[0][k];
                        xg[1*ng+i*ng2d+j*ng1d[0]+k] = x1d[1][j];
                        xg[2*ng+i*ng2d+j*ng1d[0]+k] = x1d[2][i];
                        wg[i*ng2d+j*ng1d[0]+k] = w1d[2][i]*w1d[1][j]*w1d[0][k];
                    }
            break;
        default:
            error("Dimension not implemented.\n");
    }        
}

void gaussQuadPrism(vector<double> &xg, vector<double> &wg, vector<Int> &pgauss)
{
    Int i, j, ng, nqtri, nq1d, dim = 3;
    vector<double> xtri, x1d, wtri, w1d;     
    
    gaussQuadTri(xtri, wtri,  pgauss[0]);   // triangle 
    gaussQuad1D(x1d, w1d, pgauss[1]);// 1D line    
    
    nqtri = wtri.size();
    nq1d  = w1d.size();
    ng = nq1d*nqtri;    
    xg.resize(ng*dim);
    wg.resize(ng);        
    for (i = 0; i < nq1d; i++) 
        for (j = 0; j < nqtri; j++) {
            xg[0*ng+i*nqtri+j] = xtri[0*nqtri+j];
            xg[1*ng+i*nqtri+j] = xtri[1*nqtri+j];
            xg[2*ng+i*nqtri+j] = x1d[i];
            wg[i*nqtri+j] = w1d[i]*wtri[j];
        }                
}

void gaussQuad(vector<double> &xg, vector<double> &wg, vector<Int> &pgauss, Int dim, Int elemtype)
{
    switch (dim) {
        case 0:            
            xg.resize(1);
            xg[0] = 0.0;
            wg.resize(1);
            wg[0] = 1.0;
            break;
        case 1:
            gaussQuadTensor(xg, wg, pgauss);
            break;
        case 2:
            if (elemtype == 0)
                gaussQuadTri(xg, wg, pgauss[0]);          // tri
            else if (elemtype == 1)
                gaussQuadTensor(xg, wg, pgauss);     // quad
            break;
        case 3:
            if (elemtype == 0)
                gaussQuadTet(xg, wg, pgauss[0]);     // tet
            else if (elemtype == 1)
                gaussQuadTensor(xg, wg, pgauss);     // hex
            else if (elemtype == 2)
                gaussQuadPrism(xg, wg, pgauss);     // prism
            break;
        default:
            error("Dimension not implemented.\n");
    }
}

void gaussQuad(vector<double> &xg, vector<double> &wg, Int pgauss, Int dim, Int elemtype)
{
    vector< Int > p(dim,0);
    for (Int i=0; i<dim; i++)
        p[i] = pgauss;    
    gaussQuad(xg, wg, p, dim, elemtype);
}

#endif
