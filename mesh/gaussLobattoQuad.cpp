#ifndef __GAUSSLOBATTOQUAD
#define __GAUSSLOBATTOQUAD

// Written by: C. Nguyen & P. Fernandez

void gaussLobattoQuad1d(vector<double> &x, vector<double> &w, Int pgauss)
{
    // x: nq
    // w: nq
    // Ref: http://mathworld.wolfram.com/LobattoQuadrature.html
    
    Int nq;
    if (pgauss % 2)     // pgauss is odd
        nq = (pgauss+3)/2;
    else                // pgauss is even
        nq = (pgauss+4)/2;
    
    x.resize(nq);
    w.resize(nq);
    switch (nq) {
        case 1:
            error("Gauss Lobatto quadrature not applicable for one quadrature point only...\n");
            break;
        case 2:
            x[0] = 0.0;
            x[1] = 1.0;
            w[0] = 0.5;
            w[1] = 0.5;
            break;
        case 3:
            x[0] = 0.0;
            x[1] = 0.5;
            x[2] = 1.0;
            w[0] = 0.5/3.0;
            w[1] = 2.0/3.0;
            w[2] = 0.5/3.0;
            break;
        case 4:
            x[0] =  0.0;
            x[1] = -0.5 * sqrt(5.0) / 5.0 + 0.5;
            x[2] =  0.5 * sqrt(5.0) / 5.0 + 0.5;
            x[3] =  1.0;
            w[0] = 1.0 / 12.0;
            w[1] = 5.0 / 12.0;
            w[2] = 5.0 / 12.0;
            w[3] = 1.0 / 12.0;
            break;
        case 5:
            x[0] =  0.0;
            x[1] = -0.5*sqrt(21.0) / 7.0 + 0.5;
            x[2] =  0.5;
            x[3] =  0.5*sqrt(21.0) / 7.0 + 0.5;
            x[4] =  1.0;
            w[0] = 1.0 / 20.0;
            w[1] = 49.0 / 180.0;
            w[2] = 16.0 / 45.0;
            w[3] = 49.0 / 180.0;
            w[4] = 1.0 / 20.0;
            break;
        case 6:
            x[0] =  0.0;
            x[1] = -0.5 * sqrt((7.0+2.0*sqrt(7.0)) / 21.0) + 0.5;
            x[2] = -0.5 * sqrt((7.0-2.0*sqrt(7.0)) / 21.0) + 0.5;
            x[3] =  0.5 * sqrt((7.0-2.0*sqrt(7.0)) / 21.0) + 0.5;
            x[4] =  0.5 * sqrt((7.0+2.0*sqrt(7.0)) / 21.0) + 0.5;
            x[5] =  1.0;
            w[0] = 1.0 / 30.0;
            w[1] = (14.0-sqrt(7.0)) / 60.0;
            w[2] = (14.0+sqrt(7.0)) / 60.0;
            w[3] = (14.0+sqrt(7.0)) / 60.0;
            w[4] = (14.0-sqrt(7.0)) / 60.0;
            w[5] = 1.0 / 30.0;
            break;
        default:
            error("Gauss-Lobatto quadrature not implemented for pgauss > 9.\n");
    }
}

void gaussLobattoQuadTensor(vector<double> &xg, vector<double> &wg, vector<Int> pgauss)
{
    // x: ng / nd
    // w: ng
    
    Int i, j, k, ng, ng2d, dim = pgauss.size();
    vector<Int> ng1d(dim,0);
    vector< vector<double> > x1d(dim, vector<double>());
    vector< vector<double> > w1d(dim, vector<double>());
            
    for (i=0; i<dim; i++) {        
        gaussLobattoQuad1d(x1d[i], w1d[i], pgauss[i]);
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
            for (i = 0; i < ng1d[1]; i++) {
                for (j = 0; j < ng1d[0]; j++) {
                    xg[0*ng+i*ng1d[0]+j] = x1d[0][j];
                    xg[1*ng+i*ng1d[0]+j] = x1d[1][i];
                    wg[i*ng1d[0]+j] = w1d[1][i]*w1d[0][j];
                }
            }        
            break;
        case 3:
            ng = ng1d[0]*ng1d[1]*ng1d[2];
            xg.resize(ng*3);
            wg.resize(ng);
            ng2d = ng1d[0]*ng1d[1];
            for (i = 0; i < ng1d[2]; i++) {
                for (j = 0; j < ng1d[1]; j++) {
                    for (k = 0; k < ng1d[0]; k++) {
                        xg[0*ng+i*ng2d+j*ng1d[0]+k] = x1d[0][k];
                        xg[1*ng+i*ng2d+j*ng1d[0]+k] = x1d[1][j];
                        xg[2*ng+i*ng2d+j*ng1d[0]+k] = x1d[2][i];
                        wg[i*ng2d+j*ng1d[0]+k] = w1d[2][i]*w1d[1][j]*w1d[0][k];
                    }
                }
            }   
            break;
        default:
            error("Dimension not implemented.\n");
    }        
}

void gaussLobattoQuad(vector<double> &x_p, vector<double> &w_p, vector<Int> &pgauss, Int dim, Int elemtype)
{
    printf("Gauss-Lobatto quadrature not validated yet.\n");
    
    switch (dim) {
        case 0:            
            x_p.resize(1);
            x_p[0] = 0.0;
            w_p.resize(1);
            w_p[0] = 1.0;
            break;
        case 1:
            gaussLobattoQuadTensor(x_p, w_p, pgauss);
            break;
        case 2:
            if (elemtype == 0)
                error("Gauss-Lobatto quadrature not implemented for tris.\n");
            else if (elemtype == 1)
                gaussLobattoQuadTensor(x_p, w_p, pgauss);     // quad
            break;
        case 3:
            if (elemtype == 0)
                error("Gauss-Lobatto quadrature not implemented for tets.\n");
            else if (elemtype == 1)
                gaussLobattoQuadTensor(x_p, w_p, pgauss);     // hex
            break;
        default:
            error("Dimension not implemented.\n");
    }
}

void gaussLobattoQuad(vector<double> &xg, vector<double> &wg, Int pgauss, Int dim, Int elemtype)
{
    vector< Int > p(dim,0);
    for (Int i=0; i<dim; i++)
        p[i] = pgauss;    
    gaussLobattoQuad(xg, wg, p, dim, elemtype);
}

#endif
