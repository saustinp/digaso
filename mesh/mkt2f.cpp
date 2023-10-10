#ifndef __MKT2F
#define __MKT2F

// Written by: C. Nguyen & P. Fernandez

vector< elementinfo > mkeleminfos(Int dim)
{
    Int nelemtype;
    switch (dim) {
        case 1: // line element            
            nelemtype = 1;
            break;            
        case 2:
            nelemtype = 2;
            break;
        case 3:
            nelemtype = 4;
            break;
        default:
            error("Only can handle dim=1, dim=2 or dim=3\n");
    }
    
    vector< elementinfo > eleminfos(nelemtype, elementinfo());
    
    for (Int i=0; i<nelemtype; i++)
        eleminfos[i] = mkeleminfo(dim, i);
    
    return eleminfos;    
}

vector< vector<Int> > mkf(vector< vector<Int> > &face, vector<Int> &t)
{
    
    Int i, j, nfe;
    nfe = face.size();
    
    vector< vector<Int> > f = face;    
    for (i=0; i<nfe; i++) 
        for (j = 0; j<f[i].size(); j++)
            f[i][j] = t[face[i][j]];    
    
    return f;
}


Int IsSameFace(vector<Int> face1, vector<Int> face2)
{
    Int same = 0;
    if (face1.size() == face2.size()) 
        same = IsSubset(face1,face2);
    
    return same;
}
 
// TO DO: the parallel version of this code
void mkt2f(vector<Int> &t2f, vector<Int> &t2t, vector<Int> &f, vector<Int> &t, vector<Int> elementtype, Int dim)
{
    
    vector< elementinfo> eleminfos = mkeleminfos(dim);    
    
    Int i, j, k, m, n, nve, nfe, nvf, ne, same;
    Int nvemax, nfemax, nvfmax, ntemax, elemtype, nelemtype;
    ne = elementtype.size();        
    nelemtype = eleminfos.size();   
    
    vector< Int> nfes(nelemtype,0), nves(nelemtype,0);    
    for (i=0; i<nelemtype; i++) {
        nfes[i] = eleminfos[i].nfe;
        nves[i] = eleminfos[i].nve;
    }        
    nfemax = *max_element(nfes.begin(),nfes.end());    
    ntemax = *max_element(elementtype.begin(),elementtype.end());    
    nvemax = round(((double) t.size())/ne);            
    
    vector< Int > ti(nvemax, 0);        
    vector< vector< vector<Int> > > faces(ne, vector< vector<Int> >());               
    for (i=0; i<ne; i++)
    {        
        elemtype = elementtype[i];       
        nve = nves[elemtype];        
        
        for (j=0; j<nve; j++)
            ti[j] = t[i*nvemax+j];        
        
        faces[i] = mkf(eleminfos[elemtype].face, ti);                
        
        nfe = nfes[elemtype];
        for(j = 0; j<nfe; j++)
            sort(faces[i][j].begin(), faces[i][j].end());
    }
            
    nvfmax = dim;
    if ((ntemax>0) & (dim==3))
        nvfmax = 4;            
    
    Int nmax = (1+nvfmax+2);
    t2t.resize(nfemax*ne,-1);
    t2f.resize(nfemax*ne,-1);
    f.resize(nmax*nfemax*ne,-1);        
    Int nf = 0;
    for (i=0; i<ne; i++) { // for each element i      
        nfe = nfes[elementtype[i]]; 
        nve = nves[elementtype[i]]; 
        for (k=0; k<nfe; k++) { // for each face of element i          
            if ((t2t[i*nfemax+k] > i) || (t2t[i*nfemax+k] < 0)) {
                same = 0;
                for (j=i+1; j<ne; j++) { // for each element j      
                    nfe = nfes[elementtype[j]];                        
                    for (n=0; n<nfe; n++) // for each face of element j                                      
                        if (IsSameFace(faces[i][k],faces[j][n])) {
                            same = 1;
                            break;                
                        }
                    if (same==1)
                        break;
                }
                if (same==1) {// interior face
                    t2t[i*nfemax+k] = j;
                    t2t[j*nfemax+n] = i;                    
                    t2f[i*nfemax+k] = nf;
                    t2f[j*nfemax+n] = nf;                    
                    f[nf*nmax+nmax-2] = i;
                    f[nf*nmax+nmax-1] = j;
                }
                else {// boundary face
                    t2f[i*nfemax+k] = nf;    
                    f[nf*nmax+nmax-2] = i;
                }
                f[nf*nmax] = (Int) faces[i][k].size();
                for (m = 0; m<f[nf*nmax]; m++)                     
                    f[nf*nmax+m+1] = t[i*nve+eleminfos[elementtype[i]].face[k][m]];                    
                    //f[nf*nmax+m+1] = faces[i][k][m];                                                
                nf = nf + 1;
            }
        }            
    }    
    f.resize(nmax*nf);
    
    Int reorderFaces=1;    
    if (reorderFaces == 1) {        // Reorder faces - First interior then boundary
        Int numBouFaces = 0, numNotBouFaces = 0;
        vector<Int> bouFaces(nf,-1);
        vector<Int> notBouFaces(nf,-1);
        for (i = 0; i < nf; i++) {
            if (f[(i+1)*nmax-1] == -1) {
                bouFaces[numBouFaces] = i;
                numBouFaces++;
            }
            else {
                notBouFaces[numNotBouFaces] = i;
                numNotBouFaces++;
            }
        }

        vector<Int> f_tmp = f;
        vector<Int> faceMapping(nf,-1);
        for (i = 0; i < numNotBouFaces; i++)
            for (j = 0; j < nmax; j++) {
                f_tmp[i*nmax+j] = f[notBouFaces[i]*nmax+j];
                faceMapping[notBouFaces[i]] = i;
            }
        for (i = numNotBouFaces; i < nf; i++)
            for (j = 0; j < nmax; j++) {                
                f_tmp[i*nmax+j] = f[bouFaces[i]*nmax+j];
                faceMapping[bouFaces[i]] = i;
            }
        f = f_tmp;

        // Modify t2f accordingly based on new face ordering:
        vector<Int> t2f_tmp(ne*nfe,-1);
        for (i = 0; i < ne; i++)
            for (j = 0; j < nfe; j++)
                t2f_tmp[j*ne+i] = faceMapping[t2f[j*ne+i]];        
        t2f = t2f_tmp;
    }        
}

// void localbasisface(vector<double> &philocfc, vector<double> &plocfc, Int dim, Int facetype) 
// {
//     Int npf = round(((double) plocfc.size())/dim);
//     philocfc.resize(npf*dim);
//     
//     switch (dim) {
//         case 1: // line element            
//             philocfc[0] = 1;
//             break;            
//         case 2:
//             for (i=0; i<npf; i++) {            
//                 philocfc[i] = 1 - plocfc[i];
//                 philocfc[npf+i] = plocfc[i];            
//             }    
//             break;
//         case 3:
//             if (facetype==0)
//                 for (i=0; i<npf; i++) {            
//                     philocfc[i] = 1 - plocfc[i] - plocfc[npf+i];
//                     philocfc[npf+i] = plocfc[i];
//                     philocfc[2*npf+i] = plocfc[npe+i];
//                 }    
//             else
//                 for (i=0; i<npf; i++) {  
//                     philocfc[i] = (1-plocfc[i])*(1-plocfc[npf+i]);
//                     philocfc[npf+i] = plocfc[i]*(1-plocfc[npf+i]);
//                     philocfc[2*npf+i] = plocfc[i]*plocfc[npf+i];
//                     philocfc[3*npf+i] = (1-plocfc[i])*plocfc[npf+i];
//                 }                                    
//             break;
//         default:
//             error("Only can handle dim=1, dim=2 or dim=3\n");
//     }            
// }


// % Allocate nodes
// dgnodes=zeros(npl,nd,nt);
// for dim=1:nd
//   for node=1:npv
//     dp=philocal(:,node)*p(t(:,node),dim)';
//     dgnodes(:,dim,:)=dgnodes(:,dim,:)+permute(dp,[1,3,2]);
//   end
// end

// void mkdimes(vector< Int> nfe, vector< Int> nve, Int dim)
// {
//     Int nelemtype;
//     switch (dim) {
//         case 1: // line element            
//             nelemtype = 1;
//             break;            
//         case 2:
//             nelemtype = 2;
//             break;
//         case 3:
//             nelemtype = 4;
//             break;
//         default:
//             error("Only can handle dim=1, dim=2 or dim=3\n");
//     }
//     
//     vector< Int > dime;
//     nfe.resize(nelemtype);
//     nve.resize(nelemtype);    
//     for (i=0; i<nelemtype; i++) {
//         dime = mkdime(dim, i);
//         nfe[i] = dime[0];
//         nve[i] = dime[1];
//     }        
// }
// 
// vector< vector<Int> > mkface(Int dim, Int elemtype)
// {    
//     
//     vector< Int > dime = getdime(dim, elemtype);
//     Int nfe = dime[0];
//             
//     vector< vector<Int> > face(nfe, vector<Int>());
//     
//     switch (dim) {
//         case 1: // line element
//             face[0].resize(1);
//             face[1].resize(1);
//             face[0][0] = 0;
//             face[1][0] = 1;
//             break;
//         case 2:
//             if (elemtype==0) { // triangular
//                 face[0].resize(2);
//                 face[0][0] = 1;
//                 face[0][1] = 2;
//                 
//                 face[1].resize(2);
//                 face[1][0] = 2;
//                 face[1][1] = 0;
//                 
//                 face[2].resize(2);
//                 face[2][0] = 0;
//                 face[2][1] = 1;     
//             }
//             else if (elemtype==1){ // quadrilateral
//                 //face = {0,1,2,3,1,2,3,0}; // [[1,2];[2,3];[3,4];[4,1]] - 1;
//                 face[0].resize(2);
//                 face[0][0] = 0;
//                 face[0][1] = 1;
//                 
//                 face[1].resize(2);
//                 face[1][0] = 1;
//                 face[1][1] = 2;
//                 
//                 face[2].resize(2);
//                 face[2][0] = 2;
//                 face[2][1] = 3;                            
//                 
//                 face[3].resize(2);
//                 face[3][0] = 3;
//                 face[3][1] = 0;    
//             }
//             break;
//         case 3:
//             if (elemtype==0) { // tetrahedral
//                 // [[2,3,4];[1,4,3];[1,2,4];[1,3,2]] - 1
//                 face[0].resize(3);
//                 face[0][0] = 1;
//                 face[0][1] = 2;
//                 face[0][2] = 3;
//                 
//                 face[1].resize(3);
//                 face[1][0] = 0;
//                 face[1][1] = 3;
//                 face[1][2] = 2;
//                 
//                 face[2].resize(3);
//                 face[2][0] = 0;
//                 face[2][1] = 1;
//                 face[2][2] = 3;
//                 
//                 face[3].resize(3);
//                 face[3][0] = 0;
//                 face[3][1] = 2;
//                 face[3][2] = 1;                
//             }
//             else if (elemtype==1) { // hexes
//                 //face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]] - 1;                
//                 face[0].resize(4);
//                 face[0][0] = 0;
//                 face[0][1] = 3;
//                 face[0][2] = 2;
//                 face[0][3] = 1;
//                 
//                 face[1].resize(4);
//                 face[1][0] = 4;
//                 face[1][1] = 5;
//                 face[1][2] = 6;
//                 face[1][3] = 7;
//                 
//                 face[2].resize(4);
//                 face[2][0] = 0;
//                 face[2][1] = 1;
//                 face[2][2] = 5;
//                 face[2][3] = 4;
//                 
//                 face[3].resize(4);
//                 face[3][0] = 2;
//                 face[3][1] = 3;
//                 face[3][2] = 7;                
//                 face[3][3] = 6;                
//                 
//                 face[4].resize(4);
//                 face[4][0] = 1;
//                 face[4][1] = 2;
//                 face[4][2] = 6;                
//                 face[4][3] = 5;         
//                 
//                 face[5].resize(4);
//                 face[5][0] = 3;
//                 face[5][1] = 0;
//                 face[5][2] = 4;                
//                 face[5][3] = 7;         
//                 //face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]] - 1;                
//             }
//             else if (elemtype==2) { // prism
//                 //face=[[0,2,1];[3,4,5];[1,2,5,4];[2,0,3,5];[0,1,4,3]];                
//                 face[0].resize(3);
//                 face[0][0] = 0;
//                 face[0][1] = 2;
//                 face[0][2] = 1;                
//                 
//                 face[1].resize(3);
//                 face[1][0] = 3;
//                 face[1][1] = 4;
//                 face[1][2] = 5;                
//                 
//                 face[3].resize(4);
//                 face[3][0] = 1;
//                 face[3][1] = 2;
//                 face[3][2] = 5;                
//                 face[3][3] = 4;                
//                 
//                 face[2].resize(4);
//                 face[2][0] = 2;
//                 face[2][1] = 0;
//                 face[2][2] = 3;
//                 face[2][3] = 5;                                
//                 
//                 face[4].resize(4);
//                 face[4][0] = 0;
//                 face[4][1] = 1;
//                 face[4][2] = 4;                
//                 face[4][3] = 3;                         
//                 //face=[[0,2,1];[3,4,5];[1,2,5,4];[2,0,3,5];[0,1,4,3]];                
//             }            
//             else if (elemtype==3) { // pyramid
//                 //face=[[0,3,2,1];[0,1,4];[1,2,4];[2,3,4];[3,0,4]];                
//                 face[0].resize(4);
//                 face[0][0] = 0;
//                 face[0][1] = 3;
//                 face[0][2] = 2;
//                 face[0][3] = 1;
//                 
//                 face[1].resize(3);
//                 face[1][0] = 0;
//                 face[1][1] = 1;
//                 face[1][2] = 4;
//                 
//                 face[2].resize(3);
//                 face[2][0] = 1;
//                 face[2][1] = 2;
//                 face[2][2] = 4;
//                 
//                 face[3].resize(3);
//                 face[3][0] = 2;
//                 face[3][1] = 3;
//                 face[3][2] = 4;                
//                 
//                 face[4].resize(3);
//                 face[4][0] = 3;
//                 face[4][1] = 0;
//                 face[4][2] = 4;                
//             }
//             break;
//         default:
//             error("Only can handle dim=1, dim=2 or dim=3\n");
//     }
//     
//     return face;
// }
// 
// vector< vector<Int> > mkf(vector< vector<Int> > &face, vector<Int> &t)
// {
//     
//     Int i, j, nfe;
//     nfe = face.size();
//     
//     vector< vector<Int> > f = face;    
//     for (i=0; i<nfe; i++) 
//         for (j = 0; j<f[i].size(), j++)
//             f[i][j] = t[face[i][j]];    
//     
//     return f;
// }
// 
// 
// void mkt2f(vector<Int> *t2t, vector<Int> *t, Int ne, Int nfe)
// {
//     
//     Int i, j, nve, nfe, ne, nvemax, nfemax;
//     ne = elements.size;
//     
//     
//     vector< Int> nfes, nves;    
//     mkdimes(nfes, nves, dim);
//     nfemax = *max_element(nfes.begin(),nfes.end());
//     nvemax = *max_element(nves.begin(),nves.end());
//         
//     vector< vector< vector<Int> > > faces = mkfaces(dim);
//                 
//     vector< Int > ti(nvemax, 0);
//     vector< vector< vector<Int> > > f(ne, vector< vector<Int> >());       
//     
//     //elemtypes = uniqueiarray(elements);            
//     //nelemtype = elemtypes.size();
//         
//     for (i=0; i<ne; i++)
//     {        
//         elemtype = elements[i];       
//         nve = nves[elemtype];
//         
//         for (j=0; j<nve; j++)
//             ti[j] = t[j*ne + i];        
//         f[i] = mkf(faces[elemtype], ti);                
//     }
//                 
//     //vector< vector<Int> > face;
//     //face = mkface(dim, elemtype);            
// }
// 
// void mkt2f(vector<Int> *f_p, vector<Int> *t2f_p, vector<Int> *t2t_p, vector<Int> *t_p, Int dim, Int elemtype, Int ne)
// {
//     // Note:
//     //  - The face count starts at 0 (not at 1)
//     //  - All entries in t2f are non-negative, even for boundary faces
//     
//     Int i, j, k, l, nfe, nvf, nve;
//     Int reorderFaces = 1;       // 0: No reorder. 1: Boundary first, then interior.
//     
//     Int *t2t = &t2t_p[0][0];
//     Int *t = &t_p[0][0];
//     
//     Int face[];
//     vector<Int> f_tmp;
//     vector<Int> t2f_tmp;
//     vector<Int> faceMapping;
//     vector<Int> bouFaces;
//     vector<Int> notBouFaces;
//     vector<Int> a;
//     vector<Int> b;
//     vector<Int> b_sorted;
//     vector<Int> c;
//     
//     switch (dim) {
//         case 1:
//             nfe = 2;
//             nvf = 1;
//             face = {0,1};
//             break;
//         case 2:
//             if (elemtype==0)
//                 nfe = dim+1;
//                 nvf = dim;
//                 face = {1,2,0,2,0,1};      // [[2,3];[3,1];[1,2]] - 1;
//             else if (elemtype==1)
//                 nfe = 2*dim;
//                 nvf = 2*(dim-1);
//                 face = {0,1,2,3,1,2,3,0}; // [[1,2];[2,3];[3,4];[4,1]] - 1;
//             break;
//         case 3:
//             if (elemtype==0) {
//                 nfe = dim+1;
//                 nvf = dim;
//                 face = {1,0,0,0,2,3,1,2,3,2,3,1};        // [[2,3,4];[1,4,3];[1,2,4];[1,3,2]] - 1
//             }
//             else if (elemtype==1) {
//                 nfe = 2*dim;
//                 nvf = 2*(dim-1);
//                 face = {0,4,0,2,1,3,3,5,1,3,2,0,2,6,5,7,6,4,1,7,4,6,5,7};  // [[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]] - 1;
//             }
//             break;
//         default:
//             error("Only can handle dim=1, dim=2 or dim=3\n");
//     }
//     
//     // Compute number of faces:
//     Int nf = 0;
//     for (i=0; i < ne; i++)
//         for j=0; j < nfe; j++)
//             if (t2t[j*ne+i] > i || t2t[j*ne+i] < 0) {
//                 nf++;
//             }
//     
//     // Allocate memory:
//     t2f_p[0].resize(ne*nfe);
//     Int *t2f = &t2f_p[0][0];
//     f_p[0].resize(nf*(nvf+2));
//     Int *f = &f_p[0][0];
//     for (i = 0; i < f.size(); i++)
//         f[i] = 0;
//     a.resize(nvf*nfe);
//     b.resize(nvf*nfe);
//     b_sorted.resize(nvf);
//     c.resize(nfe);
//     
//     // Compute f and t2f:
//     Int ie, jf = 0;
//     for (i=0; i < ne; i++) {
//         for j=0; j < nfe; j++) {
//             if (t2t[j*ne+i] > i || t2t[j*ne+i] < 0) {
//                 ie = t2t[j*ne+i];
//                 for (Int k = 0; k < nvf; k++)
//                     f[k*nf+jf] = t[face[k*nfe+j]*ne+i];
//                 f[nvf*nf+jf] = i;
//                 f[(nvf+1)*nf+jf] = ie;
//                 t2f[j*ne+i] = jf;
// 
//                 if (ie >= 0) {
//                     for (k = 0; k < nfe; k++)
//                         for (l = 0; l < nvf; l++)
//                             a[k*nvf+l] = t[face[l*nfe+k]*ne+ie];
//                     for (k = 0; k < nvf; k++)
//                         b[k] = f[k*nf+jf];
//                     sort(&b_sorted[0], &b[0], nvf);
//                     for (k = 0; k < nfe; k++)
//                         for (l = 0; l < nvf; l++)
//                             b[k*nvf+l] = b_sorted[l];
//                     for (k = 0; k < nfe; k++) {
//                         c[k] = 0.0;
//                         for (l = 0; l < nvf; l++)
//                             c[k] += abs(a[k*nvf+l] - b[k*nvf+l]);
//                         if (c[k] == 0)
//                             break;
//                     }
//                     t2f[k*ne+ie] = jf;
//                 }
//                 jf++;
//             }
//         }
//     }
//     if (jf != nf)
//         error("Error 5G7JH9M in mkt2f.cpp\n");
//     
//     if (reorderFaces == 1) {        // Reorder faces - First interior then boundary
//         Int numBouFaces = 0, numNotBouFaces = 0;
//         bouFaces.resize(nf);
//         notBouFaces.resize(nf);
//         for (i = 0; i < nf; i++) {
//             if (f[(nvf+1)*nf+i] == 0) {
//                 bouFaces[numBouFaces] = i;
//                 numBouFaces++;
//             }
//             else {
//                 notBouFaces[numNotBouFaces] = i;
//                 numNotBouFaces++;
//             }
//         }
// 
//         f_tmp.resize(nf*(nvf+2));
//         faceMapping.resize(nf);
//         for (i = 0; i < numBouFaces; i++)
//             for (j = 0; j < nvf+2; j++) {
//                 f_tmp[j*nf+i] = f[j*nf+bouFaces[i]];
//                 faceMapping[bouFaces[i]] = i;
//             }
//         for (i = numBouFaces; i < nf; i++)
//             for (j = 0; j < nvf+2; j++) {
//                 f_tmp[j*nf+i] = f[j*nf+notBouFaces[i]];
//                 faceMapping[notBouFaces[i]] = i;
//             }
//         for (i = numBouFaces; i < nf*(nvf+2); i++)
//             f[i] = f_tmp[i];
// 
//         // Modify t2f accordingly based on new face ordering:
//         t2f_tmp.resize(ne*nfe);
//         for (i = 0; i < ne; i++)
//             for (j = 0; j < nfe; j++)
//                 t2f_tmp[j*ne+i] = faceMapping[t2f[j*ne+i]];
//         for (i = 0; i < ne*nfe; i++)
//             t2f[i] = t2f_tmp[i];
//     }
// }
// 
#endif
