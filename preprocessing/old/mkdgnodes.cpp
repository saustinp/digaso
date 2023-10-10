#ifndef __MKDGNODES
#define __MKDGNODES

void localbasiselemface(vector<double> &philocvl,vector< vector<double> > &philocfc,vector<double> &plocvl,
        vector< vector<double> > &plocfc,Int dim, Int elemtype, Int nve) 
{
    Int i, j, npf, npe = round(((double) plocvl.size())/dim);

    philocvl.resize(npe*nve);    
    Int nfe = plocfc.size();
    philocfc.resize(nfe);        
    if (dim==1) {        
        for (i=0; i<npe; i++) {    
            philocvl[i] = 1 - plocvl[i];
            philocvl[npe+i] = plocvl[i];    
        }    
        philocfc[0].resize(1);
        philocfc[0][0] = 1;        
        philocfc[1] = philocfc[0];
    }
    else if ((dim==2) && (elemtype==0))  {    // tri
        for (i=0; i<npe; i++) {            
            philocvl[i] = 1 - plocvl[i] - plocvl[npe+i];
            philocvl[npe+i] = plocvl[i];
            philocvl[2*npe+i] = plocvl[npe+i];
        }
        for (i=0; i<nfe; i++) {      
            npf = plocfc[i].size();
            philocfc[i].resize(npf*2);
            for (j=0; j<npf; j++) {
                philocfc[i][j] = 1 - plocfc[i][j];
                philocfc[i][npf+j] = plocfc[i][j];
            }
        }
    }
    else if ((dim==2) && (elemtype==1)) {  // quad
        for (i=0; i<npe; i++) {  
            philocvl[i] = (1-plocvl[i])*(1-plocvl[npe+i]);
            philocvl[npe+i] = plocvl[i]*(1-plocvl[npe+i]);
            philocvl[2*npe+i] = plocvl[i]*plocvl[npe+i];
            philocvl[3*npe+i] = (1-plocvl[i])*plocvl[npe+i];
        }        
        for (i=0; i<nfe; i++) {      
            npf = plocfc[i].size();
            philocfc[i].resize(npf*2);
            for (j=0; j<npf; j++) {
                philocfc[i][j] = 1 - plocfc[i][j];
                philocfc[i][npf+j] = plocfc[i][j];
            }
        }
    }
    else if ((dim==3) && (elemtype==0)) { // tet
        for (i=0; i<npe; i++) {            
            philocvl[i] = 1 - plocvl[i] - plocvl[npe+i] - plocvl[2*npe+i];
            philocvl[npe+i] = plocvl[i];
            philocvl[2*npe+i] = plocvl[npe+i];
            philocvl[3*npe+i] = plocvl[2*npe+i];
        }    
        for (i=0; i<nfe; i++) {      
            npf = round(((double)plocfc[i].size())/2);
            philocfc[i].resize(npf*3);
            for (j=0; j<npf; j++) {
                philocfc[i][j] = 1 - plocfc[i][j] - plocfc[i][npf+j];
                philocfc[i][npf+j] = plocfc[i][j];
                philocfc[i][2*npf+j] = plocfc[i][npf+j];
            }
        }
    }
    else if ((dim==3) && (elemtype==1)) { // hex
        for (i=0; i<npe; i++) {  
            philocvl[i] = (1-plocvl[i])*(1-plocvl[npe+i])*(1-plocvl[2*npe+i]);
            philocvl[npe+i] = plocvl[i]*(1-plocvl[npe+i])*(1-plocvl[2*npe+i]);
            philocvl[2*npe+i] = plocvl[i]*plocvl[npe+i]*(1-plocvl[2*npe+i]);
            philocvl[3*npe+i] = (1-plocvl[i])*plocvl[npe+i]*(1-plocvl[2*npe+i]);
            philocvl[4*npe+i] = (1-plocvl[i])*(1-plocvl[npe+i])*(plocvl[2*npe+i]);
            philocvl[5*npe+i] = plocvl[i]*(1-plocvl[npe+i])*(plocvl[2*npe+i]);
            philocvl[6*npe+i] = plocvl[i]*plocvl[npe+i]*(plocvl[2*npe+i]);
            philocvl[7*npe+i] = (1-plocvl[i])*plocvl[npe+i]*(plocvl[2*npe+i]);
        }        
        for (i=0; i<nfe; i++) {      
            npf = round(((double)plocfc[i].size())/2);
            philocfc[i].resize(npf*4);
            for (j=0; j<npf; j++) {
                philocfc[i][j] = (1 - plocfc[i][j])*(1 - plocfc[i][npf+j]);
                philocfc[i][npf+j] = plocfc[i][j]*(1 - plocfc[i][npf+j]);
                philocfc[i][2*npf+j] = plocfc[i][j]*plocfc[i][npf+j];
                philocfc[i][3*npf+j] = (1 - plocfc[i][j])*plocfc[i][npf+j];
            }            
        }
    }
    else if ((dim==3) && (elemtype==2)) { // prism
        for (i=0; i<npe; i++) {  
            philocvl[i] = (1 - plocvl[i] - plocvl[npe+i])*(1-plocvl[2*npe+i]);
            philocvl[npe+i] = plocvl[i]*(1-plocvl[2*npe+i]);
            philocvl[2*npe+i] = plocvl[npe+i]*(1-plocvl[2*npe+i]);
            philocvl[3*npe+i] = (1 - plocvl[i] - plocvl[npe+i])*(plocvl[2*npe+i]);
            philocvl[4*npe+i] = plocvl[i]*(plocvl[2*npe+i]);
            philocvl[5*npe+i] = plocvl[npe+i]*(plocvl[2*npe+i]);
        }    
        for (i=0; i<2; i++) {      
            npf = round(((double)plocfc[i].size())/2);
            philocfc[i].resize(npf*3);
            for (j=0; j<npf; j++) {
                philocfc[i][j] = 1 - plocfc[i][j] - plocfc[i][npf+j];
                philocfc[i][npf+j] = plocfc[i][j];
                philocfc[i][2*npf+j] = plocfc[i][npf+j];
            }
        }    
        for (i=2; i<nfe; i++) {      
            npf = round(((double)plocfc[i].size())/2);
            philocfc[i].resize(npf*4);
            for (j=0; j<npf; j++) {
                philocfc[i][j] = (1 - plocfc[i][j])*(1 - plocfc[i][npf+j]);
                philocfc[i][npf+j] = plocfc[i][j]*(1 - plocfc[i][npf+j]);
                philocfc[i][2*npf+j] = plocfc[i][j]*plocfc[i][npf+j];
                philocfc[i][3*npf+j] = (1 - plocfc[i][j])*plocfc[i][npf+j];
            }            
        }
    }
    else if ((dim==3) && (elemtype==2)) { // pyramid
        for (i=0; i<npe; i++) {  
            philocvl[i] = (1-plocvl[i])*(1-plocvl[npe+i])*(1-plocvl[2*npe+i]);
            philocvl[npe+i] = plocvl[i]*(1-plocvl[npe+i])*(1-plocvl[2*npe+i]);
            philocvl[2*npe+i] = plocvl[i]*plocvl[npe+i]*(1-plocvl[2*npe+i]);
            philocvl[3*npe+i] = (1-plocvl[i])*plocvl[npe+i]*(1-plocvl[2*npe+i]);                
            philocvl[4*npe+i] = plocvl[2*npe+i];        
        }
        i = 0;
        npf = round(((double)plocfc[i].size())/2);
        philocfc[i].resize(npf*4);
        for (j=0; j<npf; j++) {
            philocfc[i][j] = (1 - plocfc[i][j])*(1 - plocfc[i][npf+j]);
            philocfc[i][npf+j] = plocfc[i][j]*(1 - plocfc[i][npf+j]);
            philocfc[i][2*npf+j] = plocfc[i][j]*plocfc[i][npf+j];
            philocfc[i][3*npf+j] = (1 - plocfc[i][j])*plocfc[i][npf+j];
        }
        for (i=1; i<nfe; i++) {      
            npf = round(((double)plocfc[i].size())/2);
            philocfc[i].resize(npf*3);
            for (j=0; j<npf; j++) {
                philocfc[i][j] = 1 - plocfc[i][j] - plocfc[i][npf+j];
                philocfc[i][npf+j] = plocfc[i][j];
                philocfc[i][2*npf+j] = plocfc[i][npf+j];
            }
        }
    }
    else
        error("Only can handle dim=1, dim=2 or dim=3\n");
        
}

void elementnodes(double *elemnodes, double *pv, vector<double> &philocal, Int dim, Int npe, Int nve, Int nvemax)
{
    Int j, k, m; 
    //Int nve = round(((double) philocal.size())/npe);            
    for (j=0; j<dim; j++)
        for (k=0; k<npe; k++) {
            elemnodes[j*npe+k] = 0;    
            for (m=0; m<nve; m++)
                elemnodes[j*npe+k] += philocal[m*npe+k]*pv[j*nvemax+m];
        }
}

void createnodes(vector<double> &dgnodes, vector<Int> &dgnodesidx, vector<double> &p, vector<Int> &t, 
        vector<Int> &elementtype, vector<elementstruct> &elements, vector< vector<double> > &philocvl, Int dim)
{    
    Int i, j, m, n, ne, ei, npei, nvei;
    Int nvemax = round(((double) t.size())/ne);            
    Int np = round(((double) p.size())/dim);            
    
    ne = elementtype.size();        
    dgnodesidx.resize(ne+1,0);
    dgnodesidx[0] = 0;
    for (i=0; i<ne; i++) {
        ei = elementtype[i];
        npei = elements[ei].npe;        
        dgnodesidx[i+1] = dgnodesidx[i] + npei;
    }    
        
    vector<double> pv(nvemax*dim);
    dgnodes.resize(dgnodesidx[ne]*dim);    
    for (i=0; i<ne; i++) {
        ei = elementtype[i];
        npei = elements[ei].npe;        
        nvei = elements[ei].nve;              
        for (m=0; m<nvei; m++)             
            for (n=0; n<dim; n++)
                pv[n*nvemax+m] = p[n*np+t[i*nvemax+m]];        
        elementnodes(&dgnodes[dgnodesidx[i]], &pv[0], philocvl[ei], dim, npei, nvei, nvemax);
    }            
}

#endif
