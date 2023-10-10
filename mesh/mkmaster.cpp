#ifndef __MKMASTER
#define __MKMASTER

#include "mkshape.cpp"
#include "quadrature.cpp"
//#include "mkmasternodes.cpp"

// Written by: C. Nguyen & P. Fernandez
//mkmaster(master, porder, &pgauss[0], &quadType[0], nd, elemtype, nodetype);
void mkmaster(masterstruct &master, elementstruct &element, vector<Int> &pgauss, vector<Int> &quadType)
{
    Int i, j, k, l, nfe;
    
    master.pgauss = pgauss;
    Int pgaussR = pgauss[0];
    Int pgaussJ = pgauss[1];
    Int pgaussQ = pgauss[2];
    master.pgaussR = pgaussR;
    master.pgaussJ = pgaussJ;
    master.pgaussQ = pgaussQ;
    master.quadTypeR = quadType[0];
    master.quadTypeJ = quadType[1];
    master.quadTypeQ = quadType[2];    
    
    //masternodes(&element.plocvl, &element.plocfc, &element.perm, &element.npe, &element.npf, porder, dim, elemtype, nodetype);
    // High-order nodes on the master element and face:
    //masternodes(&master.plocvl, &master.plocfc, &master.npv, &master.npf, master.porder, master.nd, master.elemtype, master.nodetype);
   // Int npv = master.npv;
    //Int npf = master.npf;
    Int npv = element.npe;
    Int elemtype = element.elemtype;
    Int nd = element.dim;
    Int porder = element.porder[0];    
    
    // Quadrature rule on the master element:
    quadrature(master.gpvlR, master.gwvlR, master.nqvR, master.pgaussR, master.quadTypeR, nd, elemtype);
    quadrature(master.gpvlJ, master.gwvlJ, master.nqvJ, master.pgaussJ, master.quadTypeJ, nd, elemtype);
    quadrature(master.gpvlQ, master.gwvlQ, master.nqvQ, master.pgaussQ, master.quadTypeQ, nd, elemtype);
    Int nqvR = master.nqvR;
    Int nqvJ = master.nqvJ;
    Int nqvQ = master.nqvQ;
        
    // Shape functions and derivatives on the master element:
    mkshape(master.shapvlR, element.plocvl, master.gpvlR, nqvR, elemtype, porder, nd, npv);
    master.shapvtR.resize(npv*nqvR*(nd+1));
    master.shapvgR.resize(npv*(nd+1)*nqvR);
    master.shapvgdotshapvlR.resize(npv*npv*(nd+1)*nqvR);
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvR; j++)
            for (k = 0; k < npv; k++) {
                master.shapvtR[i*npv*nqvR+k*nqvR+j] = master.shapvlR[i*nqvR*npv+j*npv+k];
                master.shapvgR[i*nqvR*npv+j*npv+k] = master.gwvlR[j] * master.shapvlR[i*nqvR*npv+j*npv+k];
            }
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvR; j++)
            for (k = 0; k < npv; k++)
                for (l = 0; l < npv; l++)
                    master.shapvgdotshapvlR[i*nqvR*npv*npv + j*npv*npv + k*npv + l] = master.shapvgR[i*nqvR*npv + j*npv + l] * master.shapvlR[0*nqvR*npv + j*npv + k];
    
    mkshape(master.shapnv, element.plocvl, element.plocvl, npv, elemtype, porder, nd, npv);
    master.shapnvt.resize(npv*npv*(nd+1));
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < npv; j++)
            for (k = 0; k < npv; k++) {
                master.shapnvt[i*npv*npv+k*npv+j] = master.shapnv[i*npv*npv+j*npv+k];
            }
    
    mkshape(master.shapvlJ, element.plocvl, master.gpvlJ, nqvJ, elemtype, porder, nd, npv);
    master.shapvtJ.resize(npv*nqvJ*(nd+1));
    master.shapvgJ.resize(npv*(nd+1)*nqvJ);
    master.shapvgdotshapvlJ.resize(npv*npv*(nd+1)*nqvJ);
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvJ; j++)
            for (k = 0; k < npv; k++) {
                master.shapvtJ[i*npv*nqvJ+k*nqvJ+j] = master.shapvlJ[i*nqvJ*npv+j*npv+k];
                master.shapvgJ[i*nqvJ*npv+j*npv+k] = master.gwvlJ[j] * master.shapvlJ[i*nqvJ*npv+j*npv+k];
            }
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvJ; j++)
            for (k = 0; k < npv; k++)
                for (l = 0; l < npv; l++)
                    master.shapvgdotshapvlJ[i*nqvJ*npv*npv + j*npv*npv + k*npv + l] = master.shapvgJ[i*nqvJ*npv + j*npv + l] * master.shapvlJ[0*nqvJ*npv + j*npv + k];
    
    mkshape(master.shapvlQ, element.plocvl, master.gpvlQ, nqvQ, elemtype, porder, nd, npv);
    master.shapvtQ.resize(npv*nqvQ*(nd+1));
    master.shapvgQ.resize(npv*(nd+1)*nqvQ);
    master.shapvgdotshapvlQ.resize(npv*npv*(nd+1)*nqvQ);
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvQ; j++)
            for (k = 0; k < npv; k++) {
                master.shapvtQ[i*npv*nqvQ+k*nqvQ+j] = master.shapvlQ[i*nqvQ*npv+j*npv+k];
                master.shapvgQ[i*nqvQ*npv+j*npv+k] = master.gwvlQ[j] * master.shapvlQ[i*nqvQ*npv+j*npv+k];
            }
    for (i = 0; i < (nd+1); i++)
        for (j = 0; j < nqvQ; j++)
            for (k = 0; k < npv; k++)
                for (l = 0; l < npv; l++)
                    master.shapvgdotshapvlQ[i*nqvQ*npv*npv + j*npv*npv + k*npv + l] = master.shapvgQ[i*nqvQ*npv + j*npv + l] * master.shapvlQ[0*nqvQ*npv + j*npv + k];
    
    // Shape functions and derivatives on the master face:    
    if (nd > 1) {
        Int nqfR, nqfQ, nqfJ, npf;
        nfe = element.plocfc.size();        
        master.gpfcR.resize(nfe);
        master.gpfcJ.resize(nfe);
        master.gpfcQ.resize(nfe);
        master.gwfcR.resize(nfe);
        master.gwfcJ.resize(nfe);
        master.gwfcQ.resize(nfe);
        master.nqfR.resize(nfe);
        master.nqfJ.resize(nfe);
        master.nqfQ.resize(nfe);
        master.shapfcR.resize(nfe);
        master.shapfcJ.resize(nfe);
        master.shapfcQ.resize(nfe);
        master.shapftR.resize(nfe);
        master.shapftJ.resize(nfe);
        master.shapftQ.resize(nfe);
        master.shapfgR.resize(nfe);
        master.shapfgJ.resize(nfe);
        master.shapfgQ.resize(nfe);
        master.shapfgdotshapfcR.resize(nfe);
        master.shapfgdotshapfcJ.resize(nfe);
        master.shapfgdotshapfcQ.resize(nfe);        
        for (Int n=0; n<nfe; n++) {
            // Quadrature rule on the master face:
            quadrature(master.gpfcR[n], master.gwfcR[n], master.nqfR[n], pgaussR, master.quadTypeR, nd-1, elemtype);
            quadrature(master.gpfcJ[n], master.gwfcJ[n], master.nqfJ[n], pgaussJ, master.quadTypeJ, nd-1, elemtype);
            quadrature(master.gpfcQ[n], master.gwfcQ[n], master.nqfQ[n], pgaussQ, master.quadTypeQ, nd-1, elemtype);
            nqfR = master.nqfR[n];
            nqfJ = master.nqfJ[n];
            nqfQ = master.nqfQ[n];
            npf = element.npf[n];
            
            mkshape(master.shapfcR[n], element.plocfc[n], master.gpfcR[n], nqfR, elemtype, porder, nd-1, npf);      
            master.shapftR[n].resize(npf*nqfR*nd);
            master.shapfgR[n].resize(npf*nd*nqfR);
            master.shapfgdotshapfcR[n].resize(npf*npf*nd*nqfR);
            for (i = 0; i < nd; i++)
                for (j = 0; j < nqfR; j++)
                    for (k = 0; k < npf; k++) {
                        master.shapftR[n][i*npf*nqfR+k*nqfR+j] = master.shapfcR[n][i*nqfR*npf+j*npf+k];
                        master.shapfgR[n][i*nqfR*npf+j*npf+k] = master.gwfcR[n][j] * master.shapfcR[n][i*nqfR*npf+j*npf+k];
                    }
            for (i = 0; i < nd; i++)
                for (j = 0; j < nqfR; j++)
                    for (k = 0; k < npf; k++)
                        for (l = 0; l < npf; l++)
                            master.shapfgdotshapfcR[n][i*nqfR*npf*npf + j*npf*npf + k*npf + l] = master.shapfgR[n][i*nqfR*npf + j*npf + l] * master.shapfcR[n][0*nqfR*npf + j*npf + k];

            mkshape(master.shapfcJ[n], element.plocfc[n], master.gpfcJ[n], nqfJ, elemtype, porder, nd-1, npf);      
            master.shapftJ[n].resize(npf*nqfJ*nd);
            master.shapfgJ[n].resize(npf*nd*nqfJ);
            master.shapfgdotshapfcJ[n].resize(npf*npf*nd*nqfJ);
            for (i = 0; i < nd; i++)
                for (j = 0; j < nqfJ; j++)
                    for (k = 0; k < npf; k++) {
                        master.shapftJ[n][i*npf*nqfJ+k*nqfJ+j] = master.shapfcJ[n][i*nqfJ*npf+j*npf+k];
                        master.shapfgJ[n][i*nqfJ*npf+j*npf+k] = master.gwfcJ[n][j] * master.shapfcJ[n][i*nqfJ*npf+j*npf+k];
                    }
            for (i = 0; i < nd; i++)
                for (j = 0; j < nqfJ; j++)
                    for (k = 0; k < npf; k++)
                        for (l = 0; l < npf; l++)
                            master.shapfgdotshapfcJ[n][i*nqfJ*npf*npf + j*npf*npf + k*npf + l] = master.shapfgJ[n][i*nqfJ*npf + j*npf + l] * master.shapfcJ[n][0*nqfJ*npf + j*npf + k];

            mkshape(master.shapfcQ[n], element.plocfc[n], master.gpfcQ[n], nqfQ, elemtype, porder, nd-1, npf);      
            master.shapftQ[n].resize(npf*nqfQ*nd);
            master.shapfgQ[n].resize(npf*nd*nqfQ);
            master.shapfgdotshapfcQ[n].resize(npf*npf*nd*nqfQ);
            for (i = 0; i < nd; i++)
                for (j = 0; j < nqfQ; j++)
                    for (k = 0; k < npf; k++) {
                        master.shapftQ[n][i*npf*nqfQ+k*nqfQ+j] = master.shapfcQ[n][i*nqfQ*npf+j*npf+k];
                        master.shapfgQ[n][i*nqfQ*npf+j*npf+k] = master.gwfcQ[n][j] * master.shapfcQ[n][i*nqfQ*npf+j*npf+k];
                    }
            for (i = 0; i < nd; i++)
                for (j = 0; j < nqfQ; j++)
                    for (k = 0; k < npf; k++)
                        for (l = 0; l < npf; l++)
                            master.shapfgdotshapfcQ[n][i*nqfQ*npf*npf + j*npf*npf + k*npf + l] = master.shapfgQ[n][i*nqfQ*npf + j*npf + l] * master.shapfcQ[n][0*nqfQ*npf + j*npf + k];
        }
    }
    else {
        nfe = 2;
        master.gpfcR.resize(nfe);
        master.gpfcJ.resize(nfe);
        master.gpfcQ.resize(nfe);
        master.gwfcR.resize(nfe);
        master.gwfcJ.resize(nfe);
        master.gwfcQ.resize(nfe);
        master.nqfR.resize(nfe);
        master.nqfJ.resize(nfe);
        master.nqfQ.resize(nfe);
        master.shapfcR.resize(nfe);
        master.shapfcJ.resize(nfe);
        master.shapfcQ.resize(nfe);
        master.shapftR.resize(nfe);
        master.shapftJ.resize(nfe);
        master.shapftQ.resize(nfe);
        master.shapfgR.resize(nfe);
        master.shapfgJ.resize(nfe);
        master.shapfgQ.resize(nfe);
        master.shapfgdotshapfcR.resize(nfe);
        master.shapfgdotshapfcJ.resize(nfe);
        master.shapfgdotshapfcQ.resize(nfe);                
        for (Int n=0; n<nfe; n++) {
            master.shapfcR[n].resize(1);
            master.shapfcR[n][0] = 1;
            master.shapftR[n].resize(1);
            master.shapftR[n][0] = 1;
            master.shapfgR[n].resize(1);
            master.shapfgR[n][0] = 1;
            master.shapfgdotshapfcR[n].resize(1);
            master.shapfgdotshapfcR[n][0] = 1;
            master.shapfcJ[n].resize(1);
            master.shapfcJ[n][0] = 1;
            master.shapftJ[n].resize(1);
            master.shapftJ[n][0] = 1;
            master.shapfgJ[n].resize(1);
            master.shapfgJ[n][0] = 1;
            master.shapfgdotshapfcJ[n].resize(1);
            master.shapfgdotshapfcJ[n][0] = 1;
            master.shapfcQ[n].resize(1);
            master.shapfcQ[n][0] = 1;
            master.shapftQ[n].resize(1);
            master.shapftQ[n][0] = 1;
            master.shapfgQ[n].resize(1);
            master.shapfgQ[n][0] = 1;
            master.shapfgdotshapfcQ[n].resize(1);
            master.shapfgdotshapfcQ[n][0] = 1;
        }
    }
}

#endif

