#ifndef __MKELCON
#define __MKELCON
        
vector<Int> xiny(vector< vector<double> > &x, vector< vector<double> > &y, Int m, Int n)
{
// Determine if each row of x is a member of y
// If row i of x is a member of y and x(i,:) = y(j,:) then in(i) = j
// Else in(i) = -1    
    Int i, j, k;
    Int dim = x.size();
    //Int m = x[0].size();
    //Int n = y[0].size();        
    vector<Int> in(m,-1);    
    double d;    
    for (i=0; i<m; i++) 
        for (j=0; j<n; j++) {
            d = (x[0][i]-y[0][j])*(x[0][i]-y[0][j]);
            for (k=1; k<dim; k++)
                d += (x[k][i]-y[k][j])*(x[k][i]-y[k][j]);
            if (d<1e-12) {
                in[i] = j;
                break;
            }
        }
            
    return in;
}                

vector<Int> find(const vector<Int> &x, Int a, Int compare)
{
    Int i, j = 0;
    Int n = x.size();        
    vector<Int> idx(n,0);        
    
    switch (compare) {
        case 0: // equal to
            for (i=0; i<n; i++) 
                if (x[i] == a) {
                    idx[j] = i;
                    j += 1;
                }    
            break;
        case 1: // greater than
            for (i=0; i<n; i++) 
                if (x[i] > a) {
                    idx[j] = i;
                    j += 1;
                }    
            break;
        case -1: // less than
            for (i=0; i<n; i++) 
                if (x[i] < a) {
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

void mkelcon(vector<Int> &elcon, vector<Int> &elconidx, vector< vector<Int> > &npfidx, vector<double> &dgnodes, 
        vector<Int> &dgnodesidx, vector<elementstruct> elements, vector<double> &p, vector< vector< vector<double> > > &philocfc, 
        vector<Int> &f, vector<Int> &t2f, vector<Int> &elementtype, vector<Int> &isEDGface, Int dim, Int nvfmax, Int npfmax, Int nfemax)
{            
    Int i, j, k, m, n, e1, e2, ei, nvf, npe, nfe, npf, if1, if2, neg;
    Int nmax = nvfmax+3;
    Int nf = round(((double) f.size())/(nmax)); 
    Int np = round(((double) p.size())/(dim)); 
    Int nedg = 0;
    Int ndof = 0;
    vector<Int> fi(nvfmax,-1); 
    vector< vector<double> > pf(dim, vector<double>(nvfmax));
    vector< vector<double> > facenodes(dim, vector<double>(npfmax));
    vector< vector<double> > dg1(dim, vector<double>(npfmax));
    vector< vector<double> > dg2(dim, vector<double>(npfmax));
    vector<Int> tof(npfmax,-1), ieg(nf*npf,0);
    vector<Int> in, jn, kn;
    vector< vector<double> > edg(dim, vector<double>(nf*npfmax));
    for (i=0; i<nf; i++) {      
        nvf = f[i*nmax]; //nvf = elements[ei].face[if1].size();
        for (j=0; nvf; j++)
            fi[j] = f[i*nmax+j+1];
        e1 = f[i*nmax+nvfmax+1];
        e2 = f[i*nmax+nvfmax+2];        
        ei = elementtype[e1];        
        npe = elements[ei].npe;
        nfe = elements[ei].nfe;
        
        // location of face i on element e1
        for (j=0; j<nfe; j++)
            if (t2f[e1*nfemax+j]==i)
                break;
        if1 = j;          
        npf = elements[ei].npf[if1];
        
        // coordinates of the face i 
        for (m=0; m<nvf; m++)             
            for (n=0; n<dim; n++)
                pf[n][m] = p[n*np+fi[m]];        

        // nodal points on face i
        for (j=0; j<dim; j++)
            for (k=0; k<npf; k++) {
                facenodes[j][k] = 0;    
                for (m=0; m<nvf; m++)
                    facenodes[j][k] += philocfc[ei][if1][m*npf+k]*pf[j][m];
            }
                
        if (isEDGface[i]==0) { // hdg face
            // dof numbering on face i
            for (k=0; k<npf; k++)
                tof[k] = ndof + k;
            ndof += npf;  // number of degrees of freedom 
        }
        else { // edg face
            if (nedg==0) {
                for (j=0; j<dim; j++)
                    for (k=0; k<npf; k++)
                        edg[j][k] = facenodes[j][k]; // edg nodes on faces
                for (k=0; k<npf; k++) {
                    tof[k] = ndof + k;  // dof numbering on face i    
                    ieg[k] = tof[k];    // dof numbering of edg nodes
                }
                ndof += npf; 
                nedg += npf;
            }
            else {
                in = xiny(facenodes, edg, npf, nedg);  // find which rows of pf are in edg 
                jn = find(in,-1,0);  // new edg nodes have in==-1       
                kn = find(in,-1,1);  // old edg nodes have in>-1                  
                neg = jn.size();         // number of new edg nodes                         
                
                for (j=0; j<dim; j++)
                    for (k=0; k<neg; k++)
                        edg[j][nedg+k] = facenodes[j][jn[k]]; // update edg with new edg nodes            
                
                for (k=0; k<neg; k++)
                    tof[jn[k]] = ndof+k; // tof for new edge nodes            
                for (k=0; k<kn.size(); k++)
                    tof[kn[k]] = ieg[in[kn[k]]]; // tof for old edge nodes                          
//                 for (k=0; k<npf; k++) 
//                     if (in[k]>-1)                        
//                         tof[k] = ieg[in[k]];                
                                
                for (k=0; k<neg; k++) 
                    ieg[nedg+k] = ndof+k; // update ieg with dof numbering of new edge nodes                 
                
                ndof += neg; // update ndof                     
                nedg += neg; // update nedg                                    
            }                        
        }
        
        // dg nodes on face i from element e1
        for (j=0; j<dim; j++)
            for (k=0; k<npf; k++)                
                dg1[j][k] = dgnodes[dgnodesidx[e1]*dim+j*npe+elements[ei].perm[if1][k]];
            
        in = xiny(dg1, facenodes, npf, npf); // match facenodes to dg1    
        // assign dof numbering of face i to elcon from element e1                        
        for (k=0; k<npf; k++)                 
            elcon[elconidx[e1]+npfidx[ei][if1]+k] = tof[in[k]];
        
        if (e2>-1) {               
            ei = elementtype[e2];
            nfe = elements[ei].nfe;
            npe = elements[ei].npe;                    
            for (j=0; j<nfe; j++)
                if (t2f[e2*nfemax+j]==i)
                    break;
            if2 = j;  // location of face i on element e2
            npf = elements[ei].npf[if2];
            // dg nodes on face i from element e2
            for (j=0; j<dim; j++)
                for (k=0; k<npf; k++)                
                    dg2[j][k] = dgnodes[dgnodesidx[e2]*dim+j*npe+elements[ei].perm[if2][k]];            
            in = xiny(dg2, facenodes, npf, npf); // match facenodes to dg1    
            // assign dof numbering of face i to elcon from element e1                        
            for (k=0; k<npf; k++)                 
                elcon[elconidx[e2]+npfidx[ei][if2]+k] = tof[in[k]];            
        }                        
    }                        
}
       
#endif

// void localbasiselem(vector<double> &philocvl,vector< vector<double> > &philocfc,vector<double> &plocvl,
//         vector< vector<double> > &plocfc,Int dim, Int elemtype, Int nve) 
// {
//     Int i, j, npf, npe = round(((double) plocvl.size())/dim);
// 
//     philocvl.resize(npe*nve);    
//     Int nfe = plocfc.size();
//     philocfc.resize(nfe);        
//     if (dim==1) {        
//         for (i=0; i<npe; i++) {    
//             philocvl[i] = 1 - plocvl[i];
//             philocvl[npe+i] = plocvl[i];    
//         }    
//         philocfc[0].resize(1);
//         philocfc[0][0] = 1;        
//         philocfc[1] = philocfc[0];
//     }
//     else if ((dim==2) && (elemtype==0))  {    // tri
//         for (i=0; i<npe; i++) {            
//             philocvl[i] = 1 - plocvl[i] - plocvl[npe+i];
//             philocvl[npe+i] = plocvl[i];
//             philocvl[2*npe+i] = plocvl[npe+i];
//         }
//         for (i=0; i<nfe; i++) {      
//             npf = plocfc[i].size();
//             philocfc[i].resize(npf*2);
//             for (j=0; j<npf; j++) {
//                 philocfc[i][j] = 1 - plocfc[i][j];
//                 philocfc[i][npf+j] = plocfc[i][j];
//             }
//         }
//     }
//     else if ((dim==2) && (elemtype==1)) {  // quad
//         for (i=0; i<npe; i++) {  
//             philocvl[i] = (1-plocvl[i])*(1-plocvl[npe+i]);
//             philocvl[npe+i] = plocvl[i]*(1-plocvl[npe+i]);
//             philocvl[2*npe+i] = plocvl[i]*plocvl[npe+i];
//             philocvl[3*npe+i] = (1-plocvl[i])*plocvl[npe+i];
//         }        
//         for (i=0; i<nfe; i++) {      
//             npf = plocfc[i].size();
//             philocfc[i].resize(npf*2);
//             for (j=0; j<npf; j++) {
//                 philocfc[i][j] = 1 - plocfc[i][j];
//                 philocfc[i][npf+j] = plocfc[i][j];
//             }
//         }
//     }
//     else if ((dim==3) && (elemtype==0)) { // tet
//         for (i=0; i<npe; i++) {            
//             philocvl[i] = 1 - plocvl[i] - plocvl[npe+i] - plocvl[2*npe+i];
//             philocvl[npe+i] = plocvl[i];
//             philocvl[2*npe+i] = plocvl[npe+i];
//             philocvl[3*npe+i] = plocvl[2*npe+i];
//         }    
//         for (i=0; i<nfe; i++) {      
//             npf = round(((double)plocfc[i].size())/2);
//             philocfc[i].resize(npf*3);
//             for (j=0; j<npf; j++) {
//                 philocfc[i][j] = 1 - plocfc[i][j] - plocfc[i][npf+j];
//                 philocfc[i][npf+j] = plocfc[i][j];
//                 philocfc[i][2*npf+j] = plocfc[i][npf+j];
//             }
//         }
//     }
//     else if ((dim==3) && (elemtype==1)) { // hex
//         for (i=0; i<npe; i++) {  
//             philocvl[i] = (1-plocvl[i])*(1-plocvl[npe+i])*(1-plocvl[2*npe+i]);
//             philocvl[npe+i] = plocvl[i]*(1-plocvl[npe+i])*(1-plocvl[2*npe+i]);
//             philocvl[2*npe+i] = plocvl[i]*plocvl[npe+i]*(1-plocvl[2*npe+i]);
//             philocvl[3*npe+i] = (1-plocvl[i])*plocvl[npe+i]*(1-plocvl[2*npe+i]);
//             philocvl[4*npe+i] = (1-plocvl[i])*(1-plocvl[npe+i])*(plocvl[2*npe+i]);
//             philocvl[5*npe+i] = plocvl[i]*(1-plocvl[npe+i])*(plocvl[2*npe+i]);
//             philocvl[6*npe+i] = plocvl[i]*plocvl[npe+i]*(plocvl[2*npe+i]);
//             philocvl[7*npe+i] = (1-plocvl[i])*plocvl[npe+i]*(plocvl[2*npe+i]);
//         }        
//         for (i=0; i<nfe; i++) {      
//             npf = round(((double)plocfc[i].size())/2);
//             philocfc[i].resize(npf*4);
//             for (j=0; j<npf; j++) {
//                 philocfc[i][j] = (1 - plocfc[i][j])*(1 - plocfc[i][npf+j]);
//                 philocfc[i][npf+j] = plocfc[i][j]*(1 - plocfc[i][npf+j]);
//                 philocfc[i][2*npf+j] = plocfc[i][j]*plocfc[i][npf+j];
//                 philocfc[i][3*npf+j] = (1 - plocfc[i][j])*plocfc[i][npf+j];
//             }            
//         }
//     }
//     else if ((dim==3) && (elemtype==2)) { // prism
//         for (i=0; i<npe; i++) {  
//             philocvl[i] = (1 - plocvl[i] - plocvl[npe+i])*(1-plocvl[2*npe+i]);
//             philocvl[npe+i] = plocvl[i]*(1-plocvl[2*npe+i]);
//             philocvl[2*npe+i] = plocvl[npe+i]*(1-plocvl[2*npe+i]);
//             philocvl[3*npe+i] = (1 - plocvl[i] - plocvl[npe+i])*(plocvl[2*npe+i]);
//             philocvl[4*npe+i] = plocvl[i]*(plocvl[2*npe+i]);
//             philocvl[5*npe+i] = plocvl[npe+i]*(plocvl[2*npe+i]);
//         }    
//         for (i=0; i<2; i++) {      
//             npf = round(((double)plocfc[i].size())/2);
//             philocfc[i].resize(npf*3);
//             for (j=0; j<npf; j++) {
//                 philocfc[i][j] = 1 - plocfc[i][j] - plocfc[i][npf+j];
//                 philocfc[i][npf+j] = plocfc[i][j];
//                 philocfc[i][2*npf+j] = plocfc[i][npf+j];
//             }
//         }    
//         for (i=2; i<nfe; i++) {      
//             npf = round(((double)plocfc[i].size())/2);
//             philocfc[i].resize(npf*4);
//             for (j=0; j<npf; j++) {
//                 philocfc[i][j] = (1 - plocfc[i][j])*(1 - plocfc[i][npf+j]);
//                 philocfc[i][npf+j] = plocfc[i][j]*(1 - plocfc[i][npf+j]);
//                 philocfc[i][2*npf+j] = plocfc[i][j]*plocfc[i][npf+j];
//                 philocfc[i][3*npf+j] = (1 - plocfc[i][j])*plocfc[i][npf+j];
//             }            
//         }
//     }
//     else if ((dim==3) && (elemtype==2)) { // pyramid
//         for (i=0; i<npe; i++) {  
//             philocvl[i] = (1-plocvl[i])*(1-plocvl[npe+i])*(1-plocvl[2*npe+i]);
//             philocvl[npe+i] = plocvl[i]*(1-plocvl[npe+i])*(1-plocvl[2*npe+i]);
//             philocvl[2*npe+i] = plocvl[i]*plocvl[npe+i]*(1-plocvl[2*npe+i]);
//             philocvl[3*npe+i] = (1-plocvl[i])*plocvl[npe+i]*(1-plocvl[2*npe+i]);                
//             philocvl[4*npe+i] = plocvl[2*npe+i];        
//         }
//         i = 0;
//         npf = round(((double)plocfc[i].size())/2);
//         philocfc[i].resize(npf*4);
//         for (j=0; j<npf; j++) {
//             philocfc[i][j] = (1 - plocfc[i][j])*(1 - plocfc[i][npf+j]);
//             philocfc[i][npf+j] = plocfc[i][j]*(1 - plocfc[i][npf+j]);
//             philocfc[i][2*npf+j] = plocfc[i][j]*plocfc[i][npf+j];
//             philocfc[i][3*npf+j] = (1 - plocfc[i][j])*plocfc[i][npf+j];
//         }
//         for (i=1; i<nfe; i++) {      
//             npf = round(((double)plocfc[i].size())/2);
//             philocfc[i].resize(npf*3);
//             for (j=0; j<npf; j++) {
//                 philocfc[i][j] = 1 - plocfc[i][j] - plocfc[i][npf+j];
//                 philocfc[i][npf+j] = plocfc[i][j];
//                 philocfc[i][2*npf+j] = plocfc[i][npf+j];
//             }
//         }
//     }
//     else
//         error("Only can handle dim=1, dim=2 or dim=3\n");
//         
// }
// 
// void elementnodes(double *elemnodes, double *pv, vector<double> &philocal, Int dim, Int npe, Int nve, Int nvemax);
// {
//     Int j, k, m; 
//     //Int nve = round(((double) philocal.size())/npe);            
//     for (j=0; j<dim; j++)
//         for (k=0; k<npe; k++) {
//             elemnodes[j*npe+k] = 0;    
//             for (m=0; m<nve; m++)
//                 elemnodes[j*npe+k] += philocal[m*npe+k]*pv[j*nvemax+m];
//         }
// }
// 
// 
// // void createnodes(vector<double> &dgnodes, vector<Int> &dgnodesidx, vector<double> &p, vector<Int> &t, 
// //         vector<Int> &elemtype, vector<Int> &porder, Int dim, Int nodetype);
// // {    
// //     Int i, j, n, m, nelem, nelemtype, e, elemi, npei, nvei;
// //     
// //     vector<elementstruct> elements = mkelemstructs(dim,porder,nodetype);
// //     nelemtype = elements.size();
// //     //Int nvemax = round(((double) t.size())/ne);            
// //     
// //     vector< vector<double> > philocvl(nelemtype, vector<double>());   
// //     vector< vector< vector<double> > > philocfc(nelemtype, vector< vector<double> >());   
// //     for (i=0; i<nelemtype; i++)
// //         localbasiselem(philocvl[i], philocfc[i], elements[i].plocvl, elements[i].plocfc, dim, elements[i].elemtype, elements[i].nve);                                    
// //         
// //     vector<Int> inde, elem;
// //     vector<double> pv;
// //     elem = uniqueiarray(elementtype);            
// //     nelem = elem.size();
// //     
// //     vector<Int> npe(nelem,0), nve(nelem,0); 
// //     for (i=0; i<nelem; i++) {
// //         elemi = elem[i];
// //         npe[i] = elements[elemi].npe;
// //         nve[i] = elements[elemi].nve;
// //     }
// //     Int npemax = *max_element(npe.begin(),npe.end());  
// //     Int nvemax = *max_element(nve.begin(),nve.end());  
// //     
// //     Int ne = elementtype.size();
// //     dgnodesidx.resize(ne+1,0);
// //     dgnodesidx[0] = 0;
// //     for (i=0; i<ne; i++) {
// //         elemi = elementtype[i];
// //         npei = elements[elemi].npe;        
// //         dgnodesidx[i+1] = dgnodesidx[i] + npei;
// //     }    
// //     dgnodes.resize(dgnodesidx[ne+1]*dim);    
// //     for (i=0; i<ne; i++) {
// //         elemi = elementtype[i];
// //         npei = elements[elemi].npe;        
// //         nvei = elements[elemi].nve;
// //         pv.resize(nvei*dim);        
// //         for (m=0; m<nvei; m++) {
// //             tim = t[i*nvemax+m];   
// //             for (n=0; n<dim; n++)
// //                 pv[n*nvei+m] = p[n*np+tim];
// //         }
// //         elementnodes(&dgnodes[i*npei*dim], &pv[0], philocvl[elemi], dim, npei, nvei);
// //     }            
// // }
// 
// void createnodes(vector<double> &dgnodes, vector<Int> &dgnodesidx, vector<double> &p, vector<Int> &t, 
//         vector<Int> &elemtype, vector<elementstruct> &elements, vector< vector<double> > &philocvl, Int dim)
// {    
//     Int i, j, ne, ei, npei, nvei;
//     Int nvemax = round(((double) t.size())/ne);            
//     Int np = round(((double) p.size())/dim);            
//     
//     ne = elementtype.size();        
//     dgnodesidx.resize(ne+1,0);
//     dgnodesidx[0] = 0;
//     for (i=0; i<ne; i++) {
//         ei = elementtype[i];
//         npei = elements[ei].npe;        
//         dgnodesidx[i+1] = dgnodesidx[i] + npei;
//     }    
//         
//     vector<double> pv(nvemax*dim);
//     dgnodes.resize(dgnodesidx[ne]*dim);    
//     for (i=0; i<ne; i++) {
//         ei = elementtype[i];
//         npei = elements[ei].npe;        
//         nvei = elements[ei].nve;              
//         for (m=0; m<nvei; m++) {
//             tim = t[i*nvemax+m];   
//             for (n=0; n<dim; n++)
//                 pv[n*nvemax+m] = p[n*np+tim];
//         }        
//         elementnodes(&dgnodes[dgnodesidx[i]], &pv[0], philocvl[ei], dim, npei, nvei, nvemax);
//     }            
// }

