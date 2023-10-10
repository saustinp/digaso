#ifndef __MKMASTERNODES
#define __MKMASTERNODES

// Written by: C. Nguyen & P. Fernandez

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

void GaussLobattoNodes(double *xi, Int porder)
{
    // Gauss-Lobatto nodes on the interval [0, 1]. Ref: http://mathworld.wolfram.com/LobattoQuadrature.html
    
    printf("Gauss-Lobatto nodes not validated yet.\n");
    
    switch (porder) {
        case 0:
            xi[0] = 0.5;
            break;
        case 1:
            xi[0] = 0.0;
            xi[1] = 1.0;
            break;
        case 2:
            xi[0] = 0.0;
            xi[1] = 0.5;
            xi[2] = 1.0;
            break;
        case 3:
            xi[0] =  0.0;
            xi[1] = -0.5 * sqrt(5.0) / 5.0 + 0.5;
            xi[2] =  0.5 * sqrt(5.0) / 5.0 + 0.5;
            xi[3] =  1.0;
            break;
        case 4:
            xi[0] =  0.0;
            xi[1] = -0.5*sqrt(21.0) / 7.0 + 0.5;
            xi[2] =  0.5;
            xi[3] =  0.5*sqrt(21.0) / 7.0 + 0.5;
            xi[4] =  1.0;
            break;
        case 5:
            xi[0] =  0.0;
            xi[1] = -0.5 * sqrt((7.0+2.0*sqrt(7.0)) / 21.0) + 0.5;
            xi[2] = -0.5 * sqrt((7.0-2.0*sqrt(7.0)) / 21.0) + 0.5;
            xi[3] =  0.5 * sqrt((7.0-2.0*sqrt(7.0)) / 21.0) + 0.5;
            xi[4] =  0.5 * sqrt((7.0+2.0*sqrt(7.0)) / 21.0) + 0.5;
            xi[5] =  1.0;
            break;
        default:
            error("Gauss-Lobatto nodes not implemented for porder > 5.\n");
    }
}

void ChebyshevNodes(double *xi, Int porder)
{
    // Chebyshev nodes on the interval [0, 1]:
    
    printf("Chebyshev nodes not validated yet.\n");
    
    Int numNodes = porder + 1;
    double PI = 3.1415926535897932;
    
    if (porder == 0)
        xi[0] = 0.5;
    else {
        for (int k = 0; k < numNodes; k++) {
            xi[k] = - cos(PI * (2.0 * (double) k + 1.0) / (2.0 * (double) numNodes));
            xi[k] = 0.5 + 0.5*xi[k];
        }
    }
}

void extendedChebyshevNodes(double *xi, Int porder)
{
    // Extended Chebyshev nodes on the interval [0, 1]:
    
    printf("Extended Chebyshev nodes not validated yet.\n");
    
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
        case 1:     // Gauss-Lobatto
            GaussLobattoNodes(xi, porder);
            break;
        case 2:     // Chebyshev
            ChebyshevNodes(xi, porder);
            break;
        case 3:     // Extended Chebyshev
            extendedChebyshevNodes(xi, porder);
            break;
        default:
            error("nodetype not implemented.\n");
    }
}

void lineNodes1d(vector<double> &plocvl, Int &npv, Int porder, Int nodetype)
{
    plocvl.resize(porder+1);    
    npv = porder + 1;
    nodes1d(&plocvl[0], porder, nodetype);    
}

void nodesTensor(vector<double> &plocvl, vector< vector<double> > &plocfc, Int &npv, vector<Int> &npf, vector<Int> &porder, Int dim, Int nodetype)
{
    
    Int i, j, k, nfe;    
    
    vector< vector<double> > ploc1d(dim, vector<double>());    
    for (i=0; i<dim; i++) {
        ploc1d[i].resize(porder[i]+1); 
        nodes1d(&ploc1d[i][0], porder[i], nodetype);
    }
        
    switch (dim) {
        case 2:
            nfe = 4;
            npf.resize(nfe);
            npv = (porder[0]+1)*(porder[1]+1);
            npf[0] = porder[0]+1;
            npf[1] = porder[1]+1;
            npf[2] = npf[0];
            npf[3] = npf[1];
            plocvl.resize(npv*dim);            
            for (i = 0; i < npf[1]; i++) 
                for (j = 0; j < npf[0]; j++) {
                    plocvl[0*npv+i*npf[0]+j] = ploc1d[0][j];
                    plocvl[1*npv+i*npf[0]+j] = ploc1d[1][i];
                }           
            plocfc.resize(nfe);
            for (i = 0; i < dim; i++) 
                plocfc[i] = ploc1d[i];                   
            plocfc[2] = plocfc[0];
            plocfc[3] = plocfc[1];
            break;
        case 3:   
            nfe = 6;
            npf.resize(nfe);
            npv = (porder[0]+1)*(porder[1]+1)*(porder[2]+1);
            npf[0] = (porder[0]+1)*(porder[1]+1);
            npf[2] = (porder[0]+1)*(porder[2]+1);
            npf[4] = (porder[1]+1)*(porder[2]+1);            
            npf[1] = npf[0];
            npf[3] = npf[2];
            npf[5] = npf[4];
            plocvl.resize(npv*dim);            
            for (i = 0; i < (porder[2]+1); i++) 
                for (j = 0; j < (porder[1]+1); j++) 
                    for (k = 0; k < (porder[0]+1); k++) {
                        plocvl[0*npv+i*(porder[0]+1)*(porder[1]+1)+j*(porder[0]+1)+k] = ploc1d[0][k];
                        plocvl[1*npv+i*(porder[0]+1)*(porder[1]+1)+j*(porder[0]+1)+k] = ploc1d[1][j];
                        plocvl[2*npv+i*(porder[0]+1)*(porder[1]+1)+j*(porder[0]+1)+k] = ploc1d[2][i];
                    }
            plocfc.resize(nfe);
            plocfc[0].resize(npf[0]*(dim-1));
            plocfc[2].resize(npf[2]*(dim-1));
            plocfc[4].resize(npf[4]*(dim-1));
            for (i = 0; i < (porder[1]+1); i++) 
                for (j = 0; j < (porder[0]+1); j++) {
                    plocfc[0][0*npf[0]+i*(porder[0]+1)+j] = ploc1d[0][j];
                    plocfc[0][1*npf[0]+i*(porder[0]+1)+j] = ploc1d[1][i];
                }            
            for (i = 0; i < (porder[2]+1); i++) 
                for (j = 0; j < (porder[0]+1); j++) {
                    plocfc[2][0*npf[1]+i*(porder[1]+1)+j] = ploc1d[0][j];
                    plocfc[2][1*npf[1]+i*(porder[1]+1)+j] = ploc1d[2][i];
                }            
            for (i = 0; i < (porder[2]+1); i++) 
                for (j = 0; j < (porder[1]+1); j++) {
                    plocfc[4][0*npf[2]+i*(porder[2]+1)+j] = ploc1d[1][j];
                    plocfc[4][1*npf[2]+i*(porder[2]+1)+j] = ploc1d[2][i];
                }            
            plocfc[1] = plocfc[0];
            plocfc[3] = plocfc[2];
            plocfc[5] = plocfc[4];            
            break;
        default:
            error("Dimension not implemented.\n");
    }
}

void trinodes2d(vector<double> &plocvl, vector<double> &plocfc, Int &npv, Int &npf, Int porder, Int nodetype)
{
    Int in, dim = 2;
    Int numNodes1d = porder + 1;
    Int numNodesTri =  ((porder + 1) * (porder + 2)) / 2;
    
    npv = numNodesTri;
    npf = numNodes1d;
    
    plocvl.resize(numNodesTri*dim);
    plocfc.resize(numNodes1d*(dim-1));
    
    if (porder == 0) {
        plocvl[0] = 1.0 / 3.0;
        plocvl[1] = 1.0 / 3.0;
        plocfc[0] = 1.0 / 2.0;
    }
    else {
        switch (nodetype) {
            case 0:     // Uniform
                in = 0;
                for (int j = 0; j < numNodes1d; j++)
                    for (int i = 0; i < numNodes1d-j; i++) {
                        plocvl[0*numNodesTri+in] = ((double) i) / ((double) numNodes1d - 1.0);
                        plocvl[1*numNodesTri+in] = ((double) j) / ((double) numNodes1d - 1.0);
                        in++;
                    }
                for (int i = 0; i < numNodes1d; i++)
                    plocfc[i] = ((double) i) / ((double) numNodes1d - 1.0);
                break;
            default:
                error("nodetype not implemented for tris.\n");
        }
    }    
}

void tetnodes3d(vector<double> &plocvl, vector<double> &plocfc, Int &npv, Int &npf, Int porder, Int nodetype)
{
    Int in, dim = 3;
    Int numNodes1d = porder + 1;
    Int numNodesTri =  ((porder + 1) * (porder + 2)) / 2;
    Int numNodesTet =  ((porder + 1) * (porder + 2) * (porder + 3)) / 6;
    
    npv = numNodesTet;
    npf = numNodesTri;
    
    plocvl.resize(numNodesTet*dim);
    plocfc.resize(numNodesTri*(dim-1));
    
    if (porder == 0) {
        plocvl[0] = 1.0 / 4.0;
        plocvl[1] = 1.0 / 4.0;
        plocvl[2] = 1.0 / 4.0;
        plocfc[0] = 1.0 / 3.0;
        plocfc[1] = 1.0 / 3.0;
    }
    else {
        switch (nodetype) {
            case 0:     // Uniform
                in = 0;
                for (int k = 0; k < numNodes1d; k++)
                    for (int j = 0; j < numNodes1d-k; j++)
                        for (int i = 0; i < numNodes1d-(j+k); i++) {
                            plocvl[0*numNodesTet+in] = ((double) i) / ((double) numNodes1d - 1.0);
                            plocvl[1*numNodesTet+in] = ((double) j) / ((double) numNodes1d - 1.0);
                            plocvl[2*numNodesTet+in] = ((double) k) / ((double) numNodes1d - 1.0);
                            in++;
                        }
                in = 0;
                for (int j = 0; j < numNodes1d; j++)
                    for (int i = 0; i < numNodes1d-j; i++) {
                        plocfc[0*numNodesTri+in] = ((double) i) / ((double) numNodes1d - 1.0);
                        plocfc[1*numNodesTri+in] = ((double) j) / ((double) numNodes1d - 1.0);
                        in++;
                    }
                break;
            default:
                error("nodetype not implemented for tets.\n");
        }
    }
}


void prismnodes3d(vector<double> &plocvl, vector< vector<double> > &plocfc, Int &npv, vector<Int> &npf, vector< Int > porder, Int nodetype)
{
    nodetype = 0;
    Int i, j, dim = 3;
    Int numNodes1d = porder[1] + 1;
    Int numNodesTri =  ((porder[0] + 1) * (porder[0] + 2)) / 2;
    Int nfe = 5;
    
    npv = numNodesTri*numNodes1d;
    npf.resize(nfe);
    npf[0] = numNodesTri;
    npf[1] = numNodesTri;
    npf[2] = numNodes1d*(porder[0]+1);    
    npf[3] = numNodes1d*(porder[0]+1);    
    npf[4] = numNodes1d*(porder[0]+1);    
        
    plocvl.resize(npv*dim);
    plocfc.resize(nfe);
    plocfc[0].resize(npf[0]);
    plocfc[2].resize(npf[2]);
    
    vector<double> ploc1dz;
    ploc1dz.resize(numNodes1d); 
    nodes1d(&ploc1dz[0], porder[1], nodetype);    
    
    vector<double> ploc1dx;
    ploc1dx.resize(porder[0]+1); 
    Int n1, n2;
    trinodes2d(plocfc[0], ploc1dx, n1, n2, porder[0], nodetype);
    plocfc[1] = plocfc[0];
    
    for (i = 0; i < (porder[1]+1); i++) 
        for (j = 0; j < (porder[0]+1); j++) {
            plocfc[2][0*npf[1]+i*(porder[0]+1)+j] = ploc1dx[j];
            plocfc[2][1*npf[1]+i*(porder[0]+1)+j] = ploc1dz[i];
        }
    plocfc[3] = plocfc[2];
    plocfc[4] = plocfc[2];
    
    for (i = 0; i < (porder[1]+1); i++) 
        for (j = 0; j < npf[0]; j++) {
            plocvl[0*npv+i*npf[0]+j] = plocfc[0][0*npf[0]+j];
            plocvl[1*npv+i*npf[0]+j] = plocfc[0][1*npf[0]+j];
            plocvl[2*npv+i*npf[0]+j] = ploc1dz[i];
        }    
}

void prismnodes3d(vector<double> &plocvl, vector< vector<double> > &plocfc, Int &npv, vector<Int> &npf, Int porder, Int nodetype)
{
    vector< Int > p(2,0);
    p[0] = porder;
    p[1] = porder;    
    prismnodes3d(plocvl, plocfc, npv, npf, p, nodetype);
}

void pyramidnodes3d(vector<double> &plocvl, vector< vector<double> > &plocfc, Int &npv, vector<Int> &npf, Int porder, Int nodetype)
{
    nodetype = 0;
    Int i,j,k,m,dim = 3;
    Int numNodes1d = porder + 1;
    Int numNodesTri =  ((porder + 1) * (porder + 2)) / 2;
    Int numNodesQuad =  (porder + 1) * (porder + 1);
    Int nfe = 5;
    
    npv = (porder+1)*(porder+2)*(2*porder+3)/6;
    npf.resize(nfe);
    npf[0] = numNodesQuad;    
    npf[1] = numNodesTri;    
    npf[2] = numNodesTri;    
    npf[3] = numNodesTri;    
    npf[4] = numNodesTri;    
        
    plocvl.resize(npv*dim);
    plocfc.resize(nfe);
    plocfc[0].resize(npf[0]);
    plocfc[1].resize(npf[1]);
    
    vector<double> ploc1d;
    ploc1d.resize(numNodes1d); 
    nodes1d(&ploc1d[0], porder, nodetype);    
        
    for (i = 0; i < (porder+1); i++) 
        for (j = 0; j < (porder+1); j++) {
            plocfc[0][0*npf[0]+i*(porder+1)+j] = ploc1d[j];
            plocfc[0][1*npf[0]+i*(porder+1)+j] = ploc1d[i];
        }
    
    //void trinodes2d(vector<double> &plocvl, vector<double> &plocfc, Int &npv, Int &npf, Int porder, Int nodetype)
    
    vector<double> ploc1dx;
    ploc1dx.resize(porder+1); 
    Int n1, n2;
    trinodes2d(plocfc[1], ploc1dx, n1, n2, porder, nodetype);
    plocfc[2] = plocfc[1];     
    plocfc[3] = plocfc[1];     
    plocfc[4] = plocfc[1];     
    
    double a, x1, x2;
    for (i = 0; i < numNodes1d; i++) {
        
        a  = ((double) i)/((double) porder);
        x1 = (1-a)*0.0 + a*0.5;
        x2 = (1-a)*1.0 + a*0.5;
        ploc1dx.resize(numNodes1d-i); 
        nodes1d(&ploc1dx[0], porder-i, nodetype);
        for (j=0; j<(numNodes1d-i); j++)
            ploc1dx[j] = x1 + (x2-x1)*ploc1dx[j];
        
        for (j = 0; j < (numNodes1d-i); j++) 
            for (k = 0; k < (numNodes1d-i); k++) {                                                    
                    plocvl[0*npv+i*m+j*(numNodes1d-i)+k] = ploc1dx[k];
                    plocvl[1*npv+i*m+j*(numNodes1d-i)+k] = ploc1dx[j];
                    plocvl[2*npv+i*m+j*(numNodes1d-i)+k] = ploc1d[i];
                }                
    }
}

vector<Int> find(const vector<double> &x, double a, Int compare)
{
    Int i, j = 0;
    Int n = x.size();    
    
    vector<Int> idx(n,0);     
    
    switch (compare) {
        case 0: // equal
            for (i=0; i<n; i++) 
                if (x[i] == a) {
                    idx[j] = i;
                    j += 1;
                }    
            break;
        case 1: // less than or equal
            for (i=0; i<n; i++) 
                if (x[i] <= a) {
                    idx[j] = i;
                    j += 1;
                }    
            break;
        case 2: // greater than or equal
            for (i=0; i<n; i++) 
                if (x[i] >= a) {
                    idx[j] = i;
                    j += 1;
                }    
            break;    
        case 3: // less than
            for (i=0; i<n; i++) 
                if (x[i] < a) {
                    idx[j] = i;
                    j += 1;
                }    
            break;    
        case 4: // greater than
            for (i=0; i<n; i++) 
                if (x[i] > a) {
                    idx[j] = i;
                    j += 1;
                }    
            break;   
        default:
            error("Comparison operator not implemented.\n");    
    }
    
    idx.erase(idx.begin()+j,idx.end());    
    
    return idx;
}

void getperm(vector< vector<Int> > &perm, vector<double> &plocvl, Int dim, Int elemtype, Int npv)
{
    Int i;    
    switch (dim) {
        case 1: // 1D
            perm.resize(2);    
            perm[0].resize(1);
            perm[1].resize(1);
            perm[0][0] = 0;
            perm[1][0] = plocvl.size();
            break;
        case 2: // 2D            
            if (elemtype==0)  {   // tri                 
                perm.resize(3);    
                vector<double> x(npv,0);
                //face=[1,2] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]+plocvl[1*npv+i]-1.0);                
                perm[0] = find(x, 0.0000001, 3);
                //face=[2,0] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]);                
                perm[1] = find(x, 0.0000001, 3);
                //face=[0,1] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]);                
                perm[2] = find(x, 0.0000001, 3);
            }
            else if (elemtype==1) { // quad
                perm.resize(4);    
                vector<double> x(npv,0);
                //face=[0,1] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]);
                perm[0] = find(x, 0.0000001, 3);
                //face=[1,2] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]-1.0);
                perm[1] = find(x, 0.0000001, 3);
                //face=[2,3] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]-1.0);
                perm[2] = find(x, 0.0000001, 3);
                //face=[3,0] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]);
                perm[3] = find(x, 0.0000001, 3);
            }
            break;
        case 3: // 3D            
            if (elemtype==0) {    // tet
                perm.resize(4);    
                vector<double> x(npv,0);
                //face=[1,2,3] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]+plocvl[1*npv+i]+plocvl[2*npv+i]-1.0);
                perm[0] = find(x, 0.0000001, 3);
                //face=[0,3,2] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]-1.0);
                perm[1] = find(x, 0.0000001, 3);
                //face=[0,1,3] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]-1.0);
                perm[2] = find(x, 0.0000001, 3);
                //face=[0,2,1] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]);
                perm[3] = find(x, 0.0000001, 3);                
            }
            else if (elemtype==1) { // hex
                perm.resize(6);    
                vector<double> x(npv,0);
                //face=[0,3,2,1] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]);
                perm[0] = find(x, 0.0000001, 3);
                //face [4,5,6,7]
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]-1.0);
                perm[1] = find(x, 0.0000001, 3);
                //face [0,1,5,4]
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]);
                perm[2] = find(x, 0.0000001, 3);
                //face [2,3,7,6]  
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]-1.0);
                perm[3] = find(x, 0.0000001, 3);
                //face [1,2,6,5]
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]-1.0);
                perm[4] = find(x, 0.0000001, 3);
                //face [3,0,4,7]
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]);
                perm[5] = find(x, 0.0000001, 3);
            }
            else if (elemtype==2) { // prism
                perm.resize(5);    
                vector<double> x(npv,0);
                //face=[0,2,1] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]);
                perm[0] = find(x, 0.0000001, 3);
                //face [3,4,5]
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]-1.0);
                perm[1] = find(x, 0.0000001, 3);
                //face=[1,2,5,4] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]+plocvl[1*npv+i]-1.0);                
                perm[2] = find(x, 0.0000001, 3);
                //face=[0,3,5,2] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]);                
                perm[3] = find(x, 0.0000001, 3);
                //face=[0,1,4,3] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]);                
                perm[4] = find(x, 0.0000001, 3);                
            }
            else if (elemtype==3) {// pyramid
                perm.resize(5);    
                vector<double> x(npv,0);
                //face=[0,3,2,1] (z=0)
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]);
                perm[0] = find(x, 0.0000001, 3);
                //face [0,1,4] (z-2*y=0)
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]-2.0*plocvl[1*npv+i]);
                perm[1] = find(x, 0.0000001, 3);
                //face=[1,2,4] (x+0.5*z=1)
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]+0.5*plocvl[2*npv+i]-1.0);                
                perm[2] = find(x, 0.0000001, 3);
                //face=[2,3,4] (y+0.5*z=1)
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]+0.5*plocvl[2*npv+i]-1.0);                
                perm[3] = find(x, 0.0000001, 3);
                //face=[3,0,4] (z-2*x=0) 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]-2.0*plocvl[0*npv+i]);       
                perm[4] = find(x, 0.0000001, 3);    
            }
            break;
        default:
            error("Dimension not implemented.\n");
    }            
}

void getperm(vector< vector<Int> > &perm, vector<double> &plocvl, Int dim, Int elemtype, Int npv)
{
    Int i;    
    switch (dim) {
        case 1: // 1D
            perm.resize(2);    
            perm[0].resize(1);
            perm[1].resize(1);
            perm[0][0] = 0;
            perm[1][0] = plocvl.size();
            break;
        case 2: // 2D            
            if (elemtype==0)  {   // tri                 
                perm.resize(3);    
                vector<double> x(npv,0);
                //face=[1,2] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]+plocvl[1*npv+i]-1.0);                
                perm[0] = find(x, 0.0000001, 3);
                //face=[2,0] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]);                
                perm[1] = find(x, 0.0000001, 3);
                //face=[0,1] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]);                
                perm[2] = find(x, 0.0000001, 3);
            }
            else if (elemtype==1) { // quad
                perm.resize(4);    
                vector<double> x(npv,0);
                //face=[0,1] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]);
                perm[0] = find(x, 0.0000001, 3);
                //face=[1,2] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]-1.0);
                perm[1] = find(x, 0.0000001, 3);
                //face=[2,3] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]-1.0);
                perm[2] = find(x, 0.0000001, 3);
                //face=[3,0] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]);
                perm[3] = find(x, 0.0000001, 3);
            }
            break;
        case 3: // 3D            
            if (elemtype==0) {    // tet
                perm.resize(4);    
                vector<double> x(npv,0);
                //face=[1,2,3] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]+plocvl[1*npv+i]+plocvl[2*npv+i]-1.0);
                perm[0] = find(x, 0.0000001, 3);
                //face=[0,3,2] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]-1.0);
                perm[1] = find(x, 0.0000001, 3);
                //face=[0,1,3] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]-1.0);
                perm[2] = find(x, 0.0000001, 3);
                //face=[0,2,1] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]);
                perm[3] = find(x, 0.0000001, 3);                
            }
            else if (elemtype==1) { // hex
                perm.resize(6);    
                vector<double> x(npv,0);
                //face=[0,3,2,1] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]);
                perm[0] = find(x, 0.0000001, 3);
                //face [4,5,6,7]
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]-1.0);
                perm[1] = find(x, 0.0000001, 3);
                //face [0,1,5,4]
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]);
                perm[2] = find(x, 0.0000001, 3);
                //face [2,3,7,6]  
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]-1.0);
                perm[3] = find(x, 0.0000001, 3);
                //face [1,2,6,5]
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]-1.0);
                perm[4] = find(x, 0.0000001, 3);
                //face [3,0,4,7]
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]);
                perm[5] = find(x, 0.0000001, 3);
            }
            else if (elemtype==2) { // prism
                perm.resize(5);    
                vector<double> x(npv,0);
                //face=[0,2,1] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]);
                perm[0] = find(x, 0.0000001, 3);
                //face [3,4,5]
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]-1.0);
                perm[1] = find(x, 0.0000001, 3);
                //face=[1,2,5,4] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]+plocvl[1*npv+i]-1.0);                
                perm[2] = find(x, 0.0000001, 3);
                //face=[0,3,5,2] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]);                
                perm[3] = find(x, 0.0000001, 3);
                //face=[0,1,4,3] 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]);                
                perm[4] = find(x, 0.0000001, 3);                
            }
            else if (elemtype==3) {// pyramid
                perm.resize(5);    
                vector<double> x(npv,0);
                //face=[0,3,2,1] (z=0)
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]);
                perm[0] = find(x, 0.0000001, 3);
                //face [0,1,4] (z-2*y=0)
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]-2.0*plocvl[1*npv+i]);
                perm[1] = find(x, 0.0000001, 3);
                //face=[1,2,4] (x+0.5*z=1)
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[0*npv+i]+0.5*plocvl[2*npv+i]-1.0);                
                perm[2] = find(x, 0.0000001, 3);
                //face=[2,3,4] (y+0.5*z=1)
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[1*npv+i]+0.5*plocvl[2*npv+i]-1.0);                
                perm[3] = find(x, 0.0000001, 3);
                //face=[3,0,4] (z-2*x=0) 
                for (i=0; i<npv; i++)
                    x[i] = abs(plocvl[2*npv+i]-2.0*plocvl[0*npv+i]);       
                perm[4] = find(x, 0.0000001, 3);    
            }
            break;
        default:
            error("Dimension not implemented.\n");
    }            
}

// void getpermnode(vector<Int> &permnode, vector<Int> &porder, Int dim, Int elemtype)
// {
//     switch (dim) {
//         case 1: // 1D
//             break;
//         case 2: // 2D            
//             if (elemtype==0) {    // tri                 
//             }
//             else if (elemtype==1) // quad                
//                 nodesTensor(plocvl, plocfc, npv, npf, porder, dim, nodetype);
//             break;
//         case 3: // 3D            
//             if (elemtype==0) {    // tet
//             }
//             else if (elemtype==1)  // hex
//                 nodesTensor(plocvl, plocfc, npv, npf, porder, dim, nodetype);            
//             else if (elemtype==2) // prism
//                 prismnodes3d(plocvl, plocfc, npv, npf, porder, nodetype);
//             else if (elemtype==3) // pyramid
//                 pyramidnodes3d(plocvl, plocfc, npv, npf, porder[0], nodetype);
//             break;
//         default:
//             error("Dimension not implemented.\n");
//     }    
// }

void masternodes(vector<double> &plocvl, vector< vector<double> > &plocfc, vector< vector<Int> > &perm, Int &npv, 
        vector<Int> &npf, vector<Int> &porder, Int dim, Int elemtype, Int nodetype)
{
    switch (dim) {
        case 1: // 1D
            npf.resize(2);
            npf[0] = 1;
            npf[1] = 1;
            plocfc.resize(2);
            plocfc[0].resize(1);
            plocfc[0][0] = 0.0;
            plocfc[1] = plocfc[0];            
            lineNodes1d(plocvl, npv, porder[0], nodetype);            
            break;
        case 2: // 2D            
            if (elemtype==0) {    // tri                 
                npf.resize(3);
                plocfc.resize(3);
                trinodes2d(plocvl, plocfc[0], npv, npf[0], porder[0], nodetype);                
                npf[1] = npf[0];
                npf[2] = npf[0];
                plocfc[1] = plocfc[0];
                plocfc[2] = plocfc[0];
            }
            else if (elemtype==1) // quad                
                nodesTensor(plocvl, plocfc, npv, npf, porder, dim, nodetype);
            break;
        case 3: // 3D            
            if (elemtype==0) {    // tet
                npf.resize(4);
                plocfc.resize(4);
                tetnodes3d(plocvl, plocfc[0], npv, npf[0], porder[0], nodetype);
                npf[1] = npf[0];
                npf[2] = npf[0];
                npf[3] = npf[0];
                plocfc[1] = plocfc[0];
                plocfc[2] = plocfc[0];
                plocfc[3] = plocfc[0];
            }
            else if (elemtype==1)  // hex
                nodesTensor(plocvl, plocfc, npv, npf, porder, dim, nodetype);            
            else if (elemtype==2) // prism
                prismnodes3d(plocvl, plocfc, npv, npf, porder, nodetype);
            else if (elemtype==3) // pyramid
                pyramidnodes3d(plocvl, plocfc, npv, npf, porder[0], nodetype);
            break;
        default:
            error("Dimension not implemented.\n");
    }
    
    getperm(perm, plocvl, dim, elemtype, npv);
}

void masternodesequalorder(vector<double> &plocvl, vector< vector<double> > &plocfc, vector< vector<Int> > &perm, Int &npv, vector< Int > &npf, Int porder, Int dim, Int elemtype, Int nodetype)
{
    vector< Int > p(dim,0);
    for (Int i=0; i<dim; i++)
        p[i] = porder;    
    masternodes(plocvl, plocfc, perm, npv, npf, p, dim, elemtype, nodetype);
}

#endif
