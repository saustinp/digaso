#ifndef __MKDMD
#define __MKDMD

// #include <vector>
// #include <algorithm>
// #include <numeric>

#include "mathutilities.cpp"

void elcon2entcon(vector<Int> &rowent2elem,vector<Int> &colent2elem,vector<Int> &rowent2ent,vector<Int> &colent2ent,
    vector<Int> &ent2ind, vector<Int> &elem2ent, vector<Int> &elem, vector<Int> &elementtype, vector<Int> &nps, Int npmax)
{
        
    vector<Int> ent; 
    // sort and remove duplications
    uniqueiarray(ent,elem2ent);    
    if (ent[0]<0)
        ent.erase(ent.begin());            
    
    Int i, j, np, ind;
    Int ndof = (Int) ent.size();
    Int ne = (Int) elem.size();
    Int entmax = (Int) ent.back();
    Int elemmax =  *max_element(elem.begin(), elem.end());
                
    //vector<Int> ent2ind(entmax);
    // entity-to-index mapping
    ent2ind.resize(entmax+1,-1);      
    for (i=0; i<ndof; i++)
        ent2ind[ent[i]] = i;    
    
    // element-to-index mapping
    vector<Int> elem2ind(elemmax+1,-1);    
    for (i=0; i<ne; i++)
        elem2ind[elem[i]] = i;

//     print1iarray(&ent[0], (Int) ent.size());
//     cout<<"----"<<endl;
//     print1iarray(&ent2ind[0], (Int) ent2ind.size());
//     cout<<"----"<<endl;
//     print1iarray(&elem[0], (Int) elem.size());
//     cout<<"----"<<endl;
//     print1iarray(&elem2ind[0], (Int) elem2ind.size());
//     error("here");
        
    // store number of neighboring elements for each entity
    //vector<Int> rowent2elem;
    rowent2elem.resize(ndof,0);        
    for (i=0; i<ndof; i++)
        rowent2elem[i] = 0;
    vector<Int> elc, k;
    elc.resize(npmax);
    for (i=0; i<ne; i++) { // for each element i
        np = nps[elementtype[i]];
        for (j=0;j<np;j++) // entities on element i  
            elc[j] = elem2ent[i*npmax+j];
        uniqueiarray(k,elc,0,np); // remove duplicate entities on element i      
        for (j=0; j<k.size(); j++) {
            ind = ent2ind[k[j]];  // get entity index
            rowent2elem[ind] += 1; // increase the number of elements by 1
        }            
    }     
    //print1iarray(&rowent2elem[0], (Int) rowent2elem.size());
    
    // cummulative sum
    partial_sum(rowent2elem.begin(), rowent2elem.end(), rowent2elem.begin());
    rowent2elem.insert(rowent2elem.begin(),0);
    
    //print1iarray(&rowent2elem[0], (Int) rowent2elem.size());
//     error("here");

    //if (rowent2elem.back()>10000)
    //    error("here");
    
    // store neighboring-element indices for each entity
    colent2elem.resize(rowent2elem.back(),0);
    for (i=0; i<rowent2elem.back(); i++)
        colent2elem[i] = 0;
    vector<Int> inc(ndof,0);
    for (i=0; i<ne; i++) {
        np = nps[elementtype[i]];
        for (j=0;j<np;j++) // entities on element i  
            elc[j] = elem2ent[i*npmax+j];
        uniqueiarray(k,elc,0,np); // remove duplicate entities on element i      
        for (j=0; j<k.size(); j++) {
            ind = ent2ind[k[j]];  // get entity index
            colent2elem[rowent2elem[ind]+inc[ind]] = elem[i];
            inc[ind] += 1; // pointer to the next element
        }                    
    }
    
//     print1iarray(&colent2elem[0], (Int) colent2elem.size());
//     error("here");
    
    // maxmum number of elements per row
    Int ni, me = 0;
    for (i=0;i<ndof;i++) {
        ni = rowent2elem[i+1] - rowent2elem[i];
        me = (ni > me) ? ni : me;
    }        
    
    // store number of neighboring entities for each entity
    rowent2ent.resize(ndof,0);    
    for (i=0; i<ndof; i++)
        rowent2ent[i] = 0;
    // store neighboring-entity indices for each entity
    colent2ent.resize(ndof*np*me,0);    
    for (i=0; i<ndof*np*me; i++)
        colent2ent[i] = 0;
    Int ei, l, a, e, m = 0;
    for (i=0; i<ndof; i++) {   
        ei = ent[i];  // entity ei
        ni = rowent2elem[i+1] - rowent2elem[i];
        elc.resize(ni*npmax);
        l = 0;
        for (j=0; j<ni; j++) // for each element neighboring the entity ei
        {
            e = colent2elem[rowent2elem[i]+j]; // element neighboring the entity ei
            e = elem2ind[e];
            for (a=0;a<np;a++) { // entities on element i  
                elc[l] = elem2ent[e*np+a];                        
                l += 1; 
            }
        }
        uniqueiarray(k,elc,0,l); // remove duplicate entities
        l = (Int) k.size();
        rowent2ent[i] = l; // number of neighboring entities        
        colent2ent[m] = ei;  // the entity ei
        // store its neighbors
        a = 1;
        for (j=0; j<l; j++)        
            if (k[j] != ei) {
                colent2ent[m+a] = k[j];
                a += 1;
            }        
        m += l;
    }
    // cummulative sum
    partial_sum(rowent2ent.begin(), rowent2ent.end(), rowent2ent.begin());
    rowent2ent.insert(rowent2ent.begin(),0);
    colent2ent.erase(colent2ent.begin()+m,colent2ent.end()); 
    
//     print1iarray(&rowent2ent[0], (Int) rowent2ent.size());
//     print1iarray(&colent2ent[0], (Int) colent2ent.size());
//     error("here");
}

void elcon2entcon(vector<Int> &rowent2elem,vector<Int> &colent2elem,vector<Int> &rowent2ent,vector<Int> &colent2ent,
       vector<Int> &elem2ent, vector<Int> &elementtype, vector<Int> &nps, Int npmax, Int ne)
{
    vector<Int> ent; 
    // sort and remove duplications
    uniqueiarray(ent,elem2ent);    
    if (ent[0]<0)
        ent.erase(ent.begin());            

    Int i, j, np, ind;
    Int ndof = (Int) ent.size();    
    Int entmax = (Int) ent.back();        
    
    // store number of neighboring elements for each entity    
    rowent2elem.resize(ndof,0);
    for (i=0; i<ndof; i++)
        rowent2elem[i] = 0;
    
    vector<Int> elc, k;
    elc.resize(npmax);
    for (i=0; i<ne; i++) { // for each element i
        np = nps[elementtype[i]];
        for (j=0;j<np;j++) // entities on element i  
            elc[j] = elem2ent[i*npmax+j];
        uniqueiarray(k,elc,0,np); // remove duplicate entities on element i      
        for (j=0; j<k.size(); j++) {
            ind = k[j];  // get entity index
            rowent2elem[ind] += 1; // increase the number of elements by 1
        }            
    }     
    // cummulative sum
    partial_sum(rowent2elem.begin(), rowent2elem.end(), rowent2elem.begin());
    rowent2elem.insert(rowent2elem.begin(),0);
    
    // store neighboring-element indices for each entity
    colent2elem.resize(rowent2elem.back(),0);
    for (i=0; i<rowent2elem.back(); i++)
        colent2elem[i] = 0;
    
    vector<Int> inc(ndof,0);
    for (i=0; i<ne; i++) {
        np = nps[elementtype[i]];
        for (j=0;j<np;j++) // entities on element i  
            elc[j] = elem2ent[i*npmax+j];
        uniqueiarray(k,elc,0,np); // remove duplicate entities on element i      
        for (j=0; j<k.size(); j++) {
            ind = k[j];  // get entity index
            colent2elem[rowent2elem[ind]+inc[ind]] = i;
            inc[ind] += 1; // pointer to the next element
        }                    
    }    
    
    Int ni, me = 0;
    for (i=0;i<ndof;i++) {
        ni = rowent2elem[i+1] - rowent2elem[i];
        me = (ni > me) ? ni : me;
    }        
    
    // store number of neighboring entities for each entity
    rowent2ent.resize(ndof,0);    
    for (i=0; i<ndof; i++)
        rowent2ent[i] = 0;
    // store neighboring-entity indices for each entity
    colent2ent.resize(ndof*npmax*me,0);    
    for (i=0; i<ndof*npmax*me; i++)
        colent2ent[i] = 0;
    Int ei, l, a, e, m = 0;
    for (i=0; i<ndof; i++) {   
        ei = ent[i];  // entity ei
        ni = rowent2elem[i+1] - rowent2elem[i];
        elc.resize(ni*npmax);
        l = 0;
        for (j=0; j<ni; j++) // for each element neighboring the entity ei
        {
            e = colent2elem[rowent2elem[i]+j]; // element neighboring the entity ei
            np = nps[elementtype[e]];
            for (a=0;a<np;a++) {// entities on element i  
                elc[l] = elem2ent[e*npmax+a];                        
                l += 1;
            }
        }
        uniqueiarray(k,elc,0,l); // remove duplicate entities
        l = (Int) k.size();
        rowent2ent[i] = l; // number of neighboring entities        
        colent2ent[m] = ei;  // the entity ei
        // store its neighbors
        a = 1;
        for (j=0; j<l; j++)        
            if (k[j] != ei) {
                colent2ent[m+a] = k[j];
                a += 1;
            }        
        m += l;
    }        
    // cummulative sum
    partial_sum(rowent2ent.begin(), rowent2ent.end(), rowent2ent.begin());
    rowent2ent.insert(rowent2ent.begin(),0);
    colent2ent.erase(colent2ent.begin()+m,colent2ent.end());         
}

void find_intfent(vector<Int> &intfent, vector<Int> &extent, vector<Int> &elemintfent, 
        vector<Int> &extelem, vector<Int> &intfelem, vector<Int> &intelem, vector<Int> &intent, 
        vector<Int> &ent2ind, vector<Int> &rowent2elem, vector<Int> &colent2elem, 
        vector<Int> &rowent2ent, vector<Int> &colent2ent)
{    
    Int i, j, k, nj, ij, ni, nimax, nm, m, n, n0, n1, nent, itisin;
    
    nent = intent.size();
    intfent.resize(nent);    
    for (i=0; i<nent; i++)
        intfent[i] = 0;    
    
    nimax = 0;
    for(i=0;i<nent;i++) { // For each entity in intent
        nj = intent[i]; // entity nj
        ij = ent2ind[nj]; // index of nj
        ni = rowent2ent[ij+1] - rowent2ent[ij];      
        if (ni>nimax)
            nimax = ni;
    }
    
    //find all interface entities in intent
    vector<Int> entnb(nimax,0);        
    m = 0;     
    for(i=0;i<nent;i++) { // For each entity in intent
        nj = intent[i]; // entity nj
        ij = ent2ind[nj]; // index of nj
        ni = rowent2ent[ij+1] - rowent2ent[ij];        
        //entnb.resize(ni); // list of neighboring entities 
        for (j=0; j<ni; j++) // for each entity neighboring the entity ei        
            entnb[j] = colent2ent[rowent2ent[ij]+j]; // entity neighboring the entity ei                                
        if (IsSubset(entnb, intent, ni, nent) == 0) { // entnb is not a subset of intent
            intfent[m] = intent[i];  
            m += 1;
        }                
    }
    intfent.resize(m);    

    //cout<<nent<<"  "<<m<<endl;    
    //print1iarray(&intent[0], (Int) intent.size());
    //print1iarray(&intfent[0], (Int) intfent.size());
    //print1iarray(&elemintfent[0], (Int) elemintfent.size());
        
    // find exterior entities
    vector<Int> ent0(m*nimax,0);
    n = 0;
    for(i=0;i<m;i++) { // For each entity in intent
        nj = intfent[i]; // entity nj
        ij = ent2ind[nj]; // index of nj
        ni = rowent2ent[ij+1] - rowent2ent[ij];
        //entnb.resize(ni); // list of neighboring entities 
        for (j=0; j<ni; j++) // for each entity neighboring the entity ei        
            entnb[j] = colent2ent[rowent2ent[ij]+j]; // entity neighboring the entity ei                                
        
        for (j=0; j<ni; j++) {
            itisin = 0;
            for (k=0; k<nent; k++)
                if (intent[k]==entnb[j]) {
                    itisin = 1;
                    break;
                }
            if (itisin == 0) {
                ent0[n] = entnb[j];
                n += 1;
                //ent0.push_back(entnb[j]);            
            }
        }        
    }
    ent0.resize(n);
    uniqueiarray(extent,ent0);  // exterior entities
    ent0.resize(0);
    
    //cout<<nent<<"  "<<n<<endl;    
    
    // find interface elements and exterior elements
    vector<Int> elem, elem0(m*nimax,0), elem1(m*nimax,0), elem2;
    nm = intelem.size();    
    n0 = 0;
    n1 = 0;
    for (i=0; i<m; i++) {
        nj = intfent[i]; // entity nj
        ij = ent2ind[nj]; // index of nj
        ni = rowent2elem[ij+1] - rowent2elem[ij];     
        entnb.resize(ni,0); // list of neighboring elements
        for (j=0; j<ni; j++) // for each element neighboring the entity ei        
            entnb[j] = colent2elem[rowent2elem[ij]+j]; // element neighboring the entity ei                
        if (i==0)
            elem = entnb;
        else // append entnb to elem
            elem.insert(elem.begin(), entnb.begin(), entnb.end());
                        
        for (j=0; j<ni; j++) {
            itisin = 0;
            for (k=0; k<nm; k++)
                if (intelem[k]==entnb[j]) {
                    itisin = 1;
                    break;
                }
            if (itisin == 0) {
                //elem0.push_back(entnb[j]);
                elem0[n0] = entnb[j];
                n0 += 1;
            }
            else {
                //elem1.push_back(entnb[j]);
                elem1[n1] = entnb[j];
                n1 += 1;
            }
        }                
    }    
    elem0.resize(n0);
    elem1.resize(n1);
    
    //cout<<nent<<"  "<<n0<<"  "<<n1<<endl;    
    
    // sort and remove duplications 
    uniqueiarray(extelem,elem0);  // exterior elements
    uniqueiarray(intfelem,elem1); // interface elements       
    uniqueiarray(elemintfent,elem); // elements connected to interface entities           
    
    // update interface elements
    elem = setunion(elem0, elem1);
    elem2 = setdifference(elemintfent, elem);
    intfelem.insert(intfelem.begin(), elem2.begin(), elem2.end());          
}

void find_extent(vector<Int> &extent, vector<Int> &intfelem, vector<Int> &extelem, 
        vector<Int> &intent, vector<Int> &elem2ent, Int np)
{
    Int i, j, m;
    
    // find all entities that are connected to interface elements and exterior elements        
    Int nb = intfelem.size();
    Int nc = extelem.size();            
    vector<Int> ent((nb+nc)*np);
    for (i=0;i<nb;i++) {
        m = intfelem[i]*np;        
        for (j=0; j<np; j++)
            ent[i*np+j] = elem2ent[m+j];            
    }
    for (i=0;i<nc;i++) {
        m = extelem[i]*np;        
        for (j=0; j<np; j++)
            ent[(nb+i)*np+j] = elem2ent[m+j];            
    }        
    ent = uniqueiarray(ent);        
    
    if (ent[0]<0)
        ent.erase(ent.begin());
    
    // exterior entities    
    extent = setdifference(ent,intent);                    
}


void make_elempart(vector<Int> &elempart, vector<Int> &elempartpts, vector<Int> &intelem, vector<Int> &intfelem, vector<Int> &extelem)
{
    // [interior elements, interface elements, exterior elements]  
    elempart = setdifference(intelem,intfelem);   
    elempartpts.resize(3);
    elempartpts[0] = elempart.size();    
    elempartpts[1] = intfelem.size();
    elempartpts[2] = extelem.size();
    elempart.insert(elempart.end(), intfelem.begin(), intfelem.end());    
    elempart.insert(elempart.end(), extelem.begin(), extelem.end());    
}

void find_elemrecv(vector<Int> &elemrecv, vector<Int> &extelem, vector<Int> &intelem,
                   vector<Int> &extintelem, vector<Int> &elem2cpu)
{        
    // store elements received from neighboring subdomains to assemble linear system
    Int i, k, ei;
    Int nc = extelem.size();
    Int nb = intelem.size();
    elemrecv.resize(nc*3);
    for (i=0; i<nc; i++) {
        ei = extelem[i];  // global entity numbering
        //k = elem2cpu[nb+i]; // neighboring subdomain k
        k = find(extintelem,ei)[0];
        k = elem2cpu[k];
        elemrecv[0*nc+i] = k;  // store neighboring subdomain k
        elemrecv[1*nc+i] = nb+i;  // store local entity numbering
        elemrecv[2*nc+i] = ei;  // store global entity numbering
    }                
    //sort_columnk(elemrecv, nc, 3, 0);    
    sort_rows(elemrecv, nc, 3);
}

void find_entrecv(vector<Int> &entrecv, vector<Int> &nbsd, vector<Int> &extent, vector<Int> &intent,
                  vector<Int> &extintent, vector<Int> &ent2cpu)
{    
    // store entities received from neighboring subdomains
    Int i, k, ei;
    Int nc = extent.size();
    Int nb = intent.size();
    entrecv.resize(nc*3);
    vector<Int> entrecv0(nc,0);
    for (i=0; i<nc; i++) {
        ei = extent[i];  // global entity numbering
        //k = ent2cpu[nb+i]; // neighboring subdomain k
        k = find(extintent,ei)[0];
        k = ent2cpu[k];
        entrecv[0*nc+i] = k;  // store neighboring subdomain k
        entrecv[1*nc+i] = nb+i;  // store local entity numbering
        entrecv[2*nc+i] = ei;  // store global entity numbering
        entrecv0[i] = k;  // store neighboring subdomain k
    }                
    //sort_columnk(entrecv, nc, 3, 0);
    sort_rows(entrecv, nc, 3);
    // list of neighboring subdomains from which info is required for matrix-vector product            
    nbsd = uniqueiarray(entrecv0);      
}

void find_vecrecv(vector<Int> &vecrecv, vector<Int> &extent, vector<Int> &intfent, vector<Int> &entrecv,
                  vector<Int> &ent2ind, vector<Int> &rowent2ent, vector<Int> &colent2ent)
{
    // find exterior entities connected to interface entities to perform matrix-vector product. 
    Int na = (Int) extent.size();    
    Int i, j, ni, nj, ij, k = 0;
    vector<Int> ent, entnb, on(na,0); 
    for (i=0; i<na; i++) {
        nj = entrecv[2*na+i]; // exterior global entity nj
        ij = ent2ind[nj]; // index of nj (local entity)
        ni = rowent2ent[ij+1] - rowent2ent[ij];
        entnb.resize(ni,0); // list of neirghboring entities 
        for (j=0; j<ni; j++) // for each entity neighboring the entity ei        
            entnb[j] = colent2ent[rowent2ent[ij]+j]; // entity neighboring the entity ei                 
        ent = setintersection(entnb, intfent);                
        // entity nj is connected to interface entities
        if (ent.size()>0) {
            on[i] = 1;     
            k += 1;
        }
    }    
    vecrecv.resize(k*3);
    j = 0;
    for (i=0; i<na; i++)
        if (on[i]==1) {
             vecrecv[0*k+j] = entrecv[0*na+i];  // store neighboring subdomain k
             vecrecv[1*k+j] = entrecv[1*na+i];  // store local entity numbering
             vecrecv[2*k+j] = entrecv[2*na+i];  // store global entity numbering
             j += 1;
        }
    //sort_columnk(vecrecv, k, 3, 0);        
    sort_rows(vecrecv, k, 3);
}

void find_matrecv(vector<Int> &matrecv, vector<Int> &extent, vector<Int> &elempart, vector<Int> &entrecv,
                  vector<Int> &ent2ind, vector<Int> &rowent2elem, vector<Int> &colent2elem)
{
    // find exterior entities connected to interface entities to perform matrix-vector product. 
    Int na = (Int) extent.size();    
    Int i, j, ni, nj, ij, k = 0;
    Int nelem = elempart.size();
    vector<Int> elemnb, on(na,0); 
    for (i=0; i<na; i++) {
        nj = entrecv[2*na+i]; // exterior global entity nj
        ij = ent2ind[nj]; // index of nj (local entity)
        ni = rowent2elem[ij+1] - rowent2elem[ij];
        elemnb.resize(ni); // list of neirghboring elements 
        for (j=0; j<ni; j++) // for each element neighboring the entity ei        
            elemnb[j] = colent2elem[rowent2elem[ij]+j]; // entity neighboring the entity ei                         
        // entity nj is connected to interface entities        
        if (IsSubset(elemnb,elempart,ni,nelem)==0) {
            on[i] = 1;     
            k += 1;
        }
    }    
    matrecv.resize(k*3);
        
    j = 0;
    for (i=0; i<na; i++)
        if (on[i]==1) {
             matrecv[0*k+j] = entrecv[0*na+i];  // store neighboring subdomain k
             matrecv[1*k+j] = entrecv[1*na+i];  // store local entity numbering
             matrecv[2*k+j] = entrecv[2*na+i];  // store global entity numbering
             j += 1;
        }
    //sort_columnk(matrecv, k, 3, 0);             
    sort_rows(matrecv, k, 3);
}

void make_elconhdg(vector<Int> &elcon, vector<Int> &t2f, vector<Int> &bcrs_rowent2elem,
         vector<Int> &bcrs_colent2elem, vector<Int> &bcrs_rowent2ent, vector<Int> &bcrs_colent2ent,
         vector<Int> &entpart, vector<Int> &elemmap, vector<Int> &elconhdg, vector<Int> &t2fhdg, vector<Int> elementtype, 
         vector<Int> &nfes, vector< vector<Int> > &npfs, Int nfemax, Int npfmax, Int ne)
{
    Int i, j, k, ent, nfe, npf;
    //Int ne = elementtype.size();       
    
    // local t2f for hdg
    t2f.resize(ne*nfemax);
    for (i=0;i<ne;i++) {        
        nfe = nfes[elementtype[i]];
        for(j=0;j<nfe;j++) {
            ent = t2fhdg[elemmap[i]*nfemax+j];
            for (k=0; k<entpart.size(); k++)
                if (entpart[k]==ent) {
                    t2f[i*nfemax+j] = k;  
                    break;
                }
        }        
    }
    
    // local elcon for hdg
    vector<Int> ii(npfmax), ind(npfmax);
    elcon.resize(npfmax*nfemax*ne,-1);    
    for(i=0;i<ne;i++) {
        nfe = nfes[elementtype[i]];
        for (j=0;j<nfe;j++) {
            npf = npfs[elementtype[i]][j];
            for (k=0;k<npf;k++) {
                ii[k] = elconhdg[elemmap[i]*nfemax*npfmax+j*npfmax+k] - t2fhdg[elemmap[i]*nfemax+j]*npf;
                ind[k] = t2f[i*nfemax+j]*npf+k;
            }
            for (k=0;k<npf;k++)
                elcon[i*nfemax*npfmax+j*npfmax+k] = ind[ii[k]];
        }   
    }
    
//     // take transpose
//     vector<Int> t2ft(ne*nfe);
//     for (i=0;i<ne;i++)
//         for(j=0;j<nfe;j++)         
//             t2ft[i*nfe+j] = t2f[j*ne+i];
    
    // local entity-to-element and entity-to-entity connectivities
    elcon2entcon(bcrs_rowent2elem,bcrs_colent2elem,bcrs_rowent2ent,bcrs_colent2ent,t2f,elementtype,nfes,nfemax,ne);         
}

void make_elconedg(vector<Int> &elcon, vector<Int> &t2f, vector<Int> &bcrs_rowent2elem,
         vector<Int> &bcrs_colent2elem, vector<Int> &bcrs_rowent2ent, vector<Int> &bcrs_colent2ent,
         vector<Int> &entpart, vector<Int> &elemmap, vector<Int> &elconedg, vector<Int> &t2fedg, 
         vector<Int> elementtype, vector<Int> &nfes, vector< vector<Int> > &npfs, Int nfemax, Int npfmax, Int ne)
{
    
    Int i, j, nmax, nfe, npf;
    //Int ne = (Int) elementtype.size();      
    //Int nn = round(((double) t2fedg.size()/nfe));
    t2f.resize(ne*nfemax,-1);
    for (i=0;i<ne;i++)         
        for(j=0;j<nfemax;j++)          
            t2f[i*nfemax+j] = t2fedg[elemmap[i]*nfemax+j];    
    
    vector<Int> facesInProcessor;
    uniqueiarray(facesInProcessor,t2f); 
    if (facesInProcessor[0]<0)
        facesInProcessor.erase(facesInProcessor.begin());
            
    
//      print2iarray(&t2f[0],ne,nfe);
//      printiarray(facesInProcessor);
    
    Int maxface = facesInProcessor.back();
    Int nface = (Int) facesInProcessor.size();
    vector<Int> facemapping(maxface+1,0);
    for (i=0;i<nface;i++)
        facemapping[facesInProcessor[i]] = i;
    for (i=0;i<ne;i++) {
        nfe = nfes[elementtype[i]];
        for(j=0;j<nfe;j++)          
            t2f[i*nfemax+j] = facemapping[t2f[i*nfemax+j]];
    }
    
    //print2iarray(&t2f[0],ne,nfe);
//     nmax = npfmax*nfemax;    
//     elcon.resize(nmax*ne);    
//     for(i=0;i<ne;i++) {
//         nfe = nfes[elementtype[i]];
//         npf = 0;
//         for (j=0; j<nfe; j++)
//             npf += npfs[elementtype[i]][j];                
//         for (j=0;j<npf;j++)             
//             elcon[i*nmax+j] = elconedg[elemmap[i]*nmax+j];      
//     }
    
    Int maxent = *max_element(entpart.begin(), entpart.end()); 
    Int nent = (Int) entpart.size();
    vector<Int> entmapping(maxent+1,0);
    for (i=0;i<nent;i++)
        entmapping[entpart[i]] = i;

    nmax = npfmax*nfemax;    
    elcon.resize(nmax*ne);    
    for(i=0;i<ne;i++) {
        nfe = nfes[elementtype[i]];
        npf = 0;
        for (j=0; j<nfe; j++)
            npf += npfs[elementtype[i]][j];                
        for (j=0;j<npf;j++)             
            elcon[i*nmax+j] = entmapping[elconedg[elemmap[i]*nmax+j]];
    }
    
    //print2iarray(&elcon[0],npf*nfe,ne);    
    // local entity-to-element and entity-to-entity connectivities    
    //elcon2entcon(bcrs_rowent2elem,bcrs_colent2elem,bcrs_rowent2ent,bcrs_colent2ent, elcon, nmax, ne);      
    
    vector<Int> nmaxs = nfes;
    nent = nfes.size();
    for(i=0;i<nent;i++) {
        nfe = nfes[i];
        nmaxs[i] = 0;
        for (j=0; j<nfe; j++)
            nmaxs[i] += npfs[i][j];                 
    }
    elcon2entcon(bcrs_rowent2elem,bcrs_colent2elem,bcrs_rowent2ent,bcrs_colent2ent,elcon,elementtype,nmaxs,nmax,ne);         
}

void domaindecomposition(dmdstruct &dmd, vector<Int> &elem2ent, vector<Int> &elconhdg, vector<Int> &t2f, 
        vector<Int> &elementtype, vector<Int> &extintelem, vector<Int> &extintent,  vector<Int> &elem2cpu, 
        vector<Int> &ent2cpu, vector<Int> nfes, vector< vector<Int> > npfs, 
        Int nintelem, Int nintent, Int nextelem, Int nextent, Int hdg, Int my_rank)
{   
    Int ne = elementtype.size();        
    Int nel = elem2ent.size();
    Int npmax = round(((double)nel)/ne);      
    Int nt2f = t2f.size();
    Int nfemax = round(((double) nt2f)/ne);
    Int net = elconhdg.size();
    Int npfmax = round(((double) net)/(ne*nfemax));
    
    dmd.my_rank = my_rank;        
    
    Int i, j;
    vector<Int> intelem, intfelem, extelem, intent, intfent, extent, elemintfent, elemtype;           
    intelem.resize(nintelem,0);    
    for (i=0; i<nintelem; i++)
        intelem[i] = extintelem[i];    
    intent.resize(nintent,0);
    for (i=0; i<nintent; i++)
        intent[i] = extintent[i];    
    extelem.resize(nextelem,0);
    for (i=0; i<nextelem; i++)
        extelem[i] = extintelem[nintelem+i];    
    extent.resize(nextent,0);
    for (i=0; i<nextent; i++)
        extent[i] = extintent[nintent+i];    
        
    vector<Int> nps = nfes; // HDG
    if (hdg != 1) { // not HDG
        Int nfe, n = nfes.size();        
        for(i=0;i<n;i++) {
            nfe = nfes[i];
            nps[i] = 0;
            for (j=0; j<nfe; j++)
                nps[i] += npfs[i][j];                 
        }
    }    
            
    //printf("here 1, CPU: %d\n", dmd.my_rank);
    MPI_Barrier(MPI_COMM_WORLD);  
            
    // global entity connectivities
    elcon2entcon(dmd.rowent2elem, dmd.colent2elem, dmd.rowent2ent, dmd.colent2ent,
                 dmd.ent2ind,elem2ent,extintelem,elementtype,nps,npmax);    
//     
//     if (my_rank==0)
//         printiarray(intent);
    //printf("here 2, CPU: %d\n", dmd.my_rank);
    MPI_Barrier(MPI_COMM_WORLD);  
    
    // find interface/exterior entities and elements 
    find_intfent(intfent, extent, elemintfent, extelem, intfelem, intelem, intent, 
                 dmd.ent2ind, dmd.rowent2elem, dmd.colent2elem, dmd.rowent2ent, dmd.colent2ent);
        
    // find interface elements    
    //intfelem = setintersection(elemintfent, intelem);          
  
    //printf("here 3, CPU: %d\n", dmd.my_rank);
    MPI_Barrier(MPI_COMM_WORLD);  
    
    // make element partition
    make_elempart(dmd.elempart, dmd.elempartpts, intelem, intfelem, extelem);             
   
    //printf("here 4, CPU: %d\n", dmd.my_rank);
    MPI_Barrier(MPI_COMM_WORLD);  
    
    // find the mapping from dmd.elempart to  extintelem
    Int nelem = dmd.elempart.size();
    Int nextintelem = extintelem.size();
    dmd.elemmap.resize(nelem,0);
    for (i=0; i<nelem; i++)
        for (j=0; j<nextintelem; j++) 
            if (extintelem[j]==dmd.elempart[i]) {
                dmd.elemmap[i] = j;
                break;
            }                
    elemtype = iarrayatindex(elementtype, dmd.elemmap);
    
    // make entity partition
    make_elempart(dmd.entpart, dmd.entpartpts, intent, intfent, extent);

    //printf("here 5, CPU: %d\n", dmd.my_rank);
    MPI_Barrier(MPI_COMM_WORLD);  
    
    // find the mapping from dmd.entpart to  extintent
    Int nent = dmd.entpart.size();
    Int nextintent = extintent.size();
    dmd.entmap.resize(nent,0);
    for (i=0; i<nent; i++)
        for (j=0; j<nextintent; j++) 
            if (extintent[j]==dmd.entpart[i]) {
                dmd.entmap[i] = j;
                break;
            }                    
    
    // get elemrecv from extelem
    find_elemrecv(dmd.elemrecv, extelem, intelem, extintelem, elem2cpu);
   
    // get entrecv from extent
    find_entrecv(dmd.entrecv, dmd.nbsd, extent, intent, extintent, ent2cpu);
   
    // get vecrecv
    find_vecrecv(dmd.vecrecv, extent, intfent, dmd.entrecv, dmd.ent2ind, dmd.rowent2ent, dmd.colent2ent);
   
    // get matrecv
    find_matrecv(dmd.matrecv, extent, dmd.elempart, dmd.entrecv, dmd.ent2ind, dmd.rowent2elem, dmd.colent2elem);   
    
    //printf("here 6, CPU: %d\n", dmd.my_rank);
    MPI_Barrier(MPI_COMM_WORLD);  
    
    if (hdg==1) {// hdg 
        make_elconhdg(dmd.elcon, dmd.t2f, dmd.bcrs_rowent2elem, dmd.bcrs_colent2elem, 
            dmd.bcrs_rowent2ent, dmd.bcrs_colent2ent, dmd.entpart, dmd.elemmap,
            elconhdg, t2f, elemtype, nfes, npfs, nfemax, npfmax, nelem);
    }       
    else
        make_elconedg(dmd.elcon, dmd.t2f, dmd.bcrs_rowent2elem, dmd.bcrs_colent2elem, 
            dmd.bcrs_rowent2ent, dmd.bcrs_colent2ent, dmd.entpart, dmd.elemmap,
            elem2ent, t2f, elemtype, nfes, npfs, nfemax, npfmax, nelem);    
    
    //printf("here 7, CPU: %d\n", dmd.my_rank);
    MPI_Barrier(MPI_COMM_WORLD);  
}


#ifdef HAVE_MPI

void find_datasend(vector<Int> &datasend,vector< vector<int> > &datarecv,vector<Int> &epart,vector<Int> &nbsd,Int my_rank)
{    
    Int nnbsd = (Int) nbsd.size();    
    Int nrecv, nrows, ne, e, n, i, j, k;   
    
    ne = 0;
    // determine how many rows are needed for datasend
    for (n=0; n<nnbsd; n++) { 
        nrecv = (Int) datarecv[n].size();
        nrows = round(((double)nrecv)/3);        
        for (i=0; i<nrows; i++)             
            if (datarecv[n][i]==my_rank) 
                ne = ne+1;        
    }                  
    
    datasend.resize(ne*3);
    j = 0;
    for (n=0; n<nnbsd; n++) { // for each neighboring cpu                
        nrecv = (Int) datarecv[n].size();
        nrows = round(((double)nrecv)/3);         
        for (i=0; i<nrows; i++) {            
            if (datarecv[n][i]==my_rank) {
                e = (Int) datarecv[n][2*nrows+i];
                k = find(epart,e)[0];
                datasend[0*ne+j] = nbsd[n];
                datasend[1*ne+j] = k;
                datasend[2*ne+j] = epart[k];
                j = j + 1;
            }
        }            
    }
    sort_rows(datasend, ne, 3);
}

vector< vector<int> > getdatarecv(vector<int> & datasend, vector<Int> & nbsd)
{
    // send datasend to neighboring cpus 
    // each cpu will receive and return datarecv
        
    Int n, neighbor;
    Int nnbsd = (Int) nbsd.size();
    
    MPI_Request * requests;
    MPI_Status * statuses;
    
    requests = (MPI_Request *) malloc( 2*nnbsd * sizeof(MPI_Request) );
    statuses = (MPI_Status *) malloc( 2*nnbsd * sizeof(MPI_Status) );
            
    /* Initialize request counter */
    Int requestCounter = 0;        
    Int nsend = (Int) datasend.size();
    
    // Send nsend
    int nsend_tmp[1];
    nsend_tmp[0] = (int) nsend;
    vector<int> nrecv_tmp(nnbsd, 0);            
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n]; 
        MPI_Isend(&nsend_tmp[0], 1, MPI_INT, (int) neighbor, 0,
               MPI_COMM_WORLD, &requests[requestCounter]);
            requestCounter += 1;   
    }   
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n]; 
        MPI_Irecv(&nrecv_tmp[n], 1, MPI_INT, (int) neighbor, 0,
               MPI_COMM_WORLD, &requests[requestCounter]);            
        requestCounter += 1;
    }
    
    MPI_Waitall(requestCounter, requests, statuses); 
    
    // Set requestCounter to zero:
    requestCounter = 0;    
    
    /* Non-blocking send */    
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];   
        if (nsend>0) {
            MPI_Isend(&datasend[0], (int) nsend, MPI_INT, (int) neighbor, 0,
                   MPI_COMM_WORLD, &requests[requestCounter]);
            requestCounter += 1;
        }        
    }            
    
    // initialize datarecv for all neighboring cpus   
    vector< vector<int> > datarecv(nnbsd, vector<int>());            
    
    /* Non-blocking receive */
    Int nrecv;
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];
        
        nrecv = (Int) nrecv_tmp[n];

        if (nrecv>0) {
            datarecv[n].resize(nrecv);
            MPI_Irecv(&datarecv[n][0], (int) nrecv, MPI_INT, (int) neighbor, 0,
                   MPI_COMM_WORLD, &requests[requestCounter]);            
            requestCounter += 1;
        }        
    }    
            
    /* Wait until all send and receive operations are completely done */
    MPI_Waitall(requestCounter, requests, statuses);        
    
    free(requests); 
    free(statuses);

    return datarecv;        
}

vector< vector<int> > getdatarecv(vector < vector<int> > & datasend, vector<Int> & nbsd)
{
    // send datasend to neighboring cpus 
    // each cpu will receive and return datarecv
    
    Int n, neighbor;
    Int nnbsd = (Int) nbsd.size();
    
    MPI_Request * requests;
    MPI_Status * statuses;
    requests = (MPI_Request *) malloc( 2*nnbsd * sizeof(MPI_Request) );
    statuses = (MPI_Status *) malloc( 2*nnbsd * sizeof(MPI_Status) );
            
    /* Initialize request counter */
    Int requestCounter = 0;        
    Int nsend;
    
    // Send nsend
    vector<int> nsend_tmp(nnbsd, 0);   
    vector<int> nrecv_tmp(nnbsd, 0);            
    for (n=0; n<nnbsd; n++) {
        nsend_tmp[n] = (int) datasend[n].size();
        neighbor = nbsd[n]; 
        MPI_Isend(&nsend_tmp[n], 1, MPI_INT, (int) neighbor, 0,
               MPI_COMM_WORLD, &requests[requestCounter]);
        requestCounter += 1;   
    }   
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n]; 
        MPI_Irecv(&nrecv_tmp[n], 1, MPI_INT, (int) neighbor, 0,
               MPI_COMM_WORLD, &requests[requestCounter]);            
        requestCounter += 1;
    }
    
    MPI_Waitall(requestCounter, requests, statuses); 
    
    // Set requestCounter to zero:
    requestCounter = 0;    
    
    /* Non-blocking send */    
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];   
        nsend = (Int) datasend[n].size();
        if (nsend>0) {
            MPI_Isend(&datasend[n][0], (int) nsend, MPI_INT, (int) neighbor, 0,
                   MPI_COMM_WORLD, &requests[requestCounter]);
            requestCounter += 1;
        }        
    }            
    
    // initialize datarecv for all neighboring cpus   
    vector< vector<int> > datarecv(nnbsd, vector<int>());            
    
    /* Non-blocking receive */
    Int nrecv;
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];
        
        nrecv = (Int) nrecv_tmp[n];

        if (nrecv>0) {
            datarecv[n].resize(nrecv);
            MPI_Irecv(&datarecv[n][0], (int) nrecv, MPI_INT, (int) neighbor, 0,
                   MPI_COMM_WORLD, &requests[requestCounter]);            
            requestCounter += 1;
        }
    }    
            
    /* Wait until all send and receive operations are completely done */
    MPI_Waitall(requestCounter, requests, statuses);        
    
    free(requests); 
    free(statuses);

    return datarecv;        
}

void getdatarecv(vector< vector<int> > &datarecv, vector<int> & datasend, vector<Int> & nbsd)
{
    // send datasend to neighboring cpus 
    // each cpu will receive and return datarecv
        
    Int n, neighbor;
    Int nnbsd = (Int) nbsd.size();
    
    MPI_Request * requests;
    MPI_Status * statuses;
    
    requests = (MPI_Request *) malloc( 2*nnbsd * sizeof(MPI_Request) );
    statuses = (MPI_Status *) malloc( 2*nnbsd * sizeof(MPI_Status) );
            
    /* Initialize request counter */
    Int requestCounter = 0;        
    Int nsend = (Int) datasend.size();
    
    // Send nsend
    int nsend_tmp[1];
    nsend_tmp[0] = (int) nsend;
    vector<int> nrecv_tmp(nnbsd, 0);            
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n]; 
        MPI_Isend(&nsend_tmp[0], 1, MPI_INT, (int) neighbor, 0,
               MPI_COMM_WORLD, &requests[requestCounter]);
            requestCounter += 1;   
    }   
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n]; 
        MPI_Irecv(&nrecv_tmp[n], 1, MPI_INT, (int) neighbor, 0,
               MPI_COMM_WORLD, &requests[requestCounter]);            
        requestCounter += 1;
    }
    
    MPI_Waitall(requestCounter, requests, statuses); 
    
    // Set requestCounter to zero:
    requestCounter = 0;    
    
    /* Non-blocking send */    
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];   
        if (nsend>0) {
            MPI_Isend(&datasend[0], (int) nsend, MPI_INT, (int) neighbor, 0,
                   MPI_COMM_WORLD, &requests[requestCounter]);
            requestCounter += 1;
        }        
    }            
    
    /* Non-blocking receive */
    Int nrecv;
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];
        
        nrecv = (Int) nrecv_tmp[n];

        if (nrecv>0) {
            datarecv[n].resize(nrecv);
            MPI_Irecv(&datarecv[n][0], (int) nrecv, MPI_INT, (int) neighbor, 0,
                   MPI_COMM_WORLD, &requests[requestCounter]);            
            requestCounter += 1;
        }        
    }    
            
    /* Wait until all send and receive operations are completely done */
    MPI_Waitall(requestCounter, requests, statuses);        
    
    free(requests); 
    free(statuses);
}

void getdatarecv(vector< vector<int> > & datarecv, vector < vector<int> > & datasend, vector<Int> & nbsd)
{
    // send datasend to neighboring cpus 
    // each cpu will receive and return datarecv
    
    Int n, neighbor;
    Int nnbsd = (Int) nbsd.size();
    
    MPI_Request * requests;
    MPI_Status * statuses;
    requests = (MPI_Request *) malloc( 2*nnbsd * sizeof(MPI_Request) );
    statuses = (MPI_Status *)  malloc( 2*nnbsd * sizeof(MPI_Status)  );
    
    /* Initialize request counter */
    Int requestCounter = 0;        
    Int nsend;
    
    // Send nsend
    vector<int> nsend_tmp(nnbsd, 0);   
    vector<int> nrecv_tmp(nnbsd, 0);            
    for (n=0; n<nnbsd; n++) {
        nsend_tmp[n] = (int) datasend[n].size();
        neighbor = nbsd[n]; 
        MPI_Isend(&nsend_tmp[n], 1, MPI_INT, (int) neighbor, 0,
               MPI_COMM_WORLD, &requests[requestCounter]);
        requestCounter += 1;   
    }   
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n]; 
        MPI_Irecv(&nrecv_tmp[n], 1, MPI_INT, (int) neighbor, 0,
               MPI_COMM_WORLD, &requests[requestCounter]);            
        requestCounter += 1;
    }
    
    MPI_Waitall(requestCounter, requests, statuses); 
    
    // Set requestCounter to zero:
    requestCounter = 0;    
    
    /* Non-blocking send */    
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];   
        nsend = (Int) datasend[n].size();
        if (nsend>0) {
            MPI_Isend(&datasend[n][0], (int) nsend, MPI_INT, (int) neighbor, 0,
                   MPI_COMM_WORLD, &requests[requestCounter]);
            requestCounter += 1;
        }        
    }            
    
    /* Non-blocking receive */
    Int nrecv;
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];
        
        nrecv = (Int) nrecv_tmp[n];
        
        if (nrecv>0) {
            datarecv[n].resize(nrecv);
            MPI_Irecv(&datarecv[n][0], (int) nrecv, MPI_INT, (int) neighbor, 0,
                   MPI_COMM_WORLD, &requests[requestCounter]);            
            requestCounter += 1;
        }        
    }    
            
    /* Wait until all send and receive operations are completely done */
    MPI_Waitall(requestCounter, requests, statuses);        
    
    free(requests); 
    free(statuses);     
}

void make_nbdmd(vector< dmdstruct > &nbdmd, dmdstruct &dmd)
{
    Int nnbsd = (Int) dmd.nbsd.size();                        
    vector< vector<int> > v(nnbsd, vector<int>());
    vector<int> send_tmp;
    
    send_tmp.resize(dmd.elemrecv.size());
    for (Int i=0; i<send_tmp.size(); i++)
        send_tmp[i] = (int) dmd.elemrecv[i];
    getdatarecv(v, send_tmp, dmd.nbsd);    // get elemrecv from neighbors and store it in v
    find_datasend(dmd.elemsend, v, dmd.elempart, dmd.nbsd, dmd.my_rank);    // make elemsend using v
    for (Int n=0; n<nnbsd; n++) {   // get elemrecv for my neighbors
//         nbdmd[n].elemrecv.resize(v[n].size());
//         for (Int i=0; i<v[n].size(); i++)
//             nbdmd[n].elemrecv[i] = (Int) v[n][i];
        v[n].clear();
    }
    
    send_tmp.resize(dmd.entrecv.size());
    for (Int i=0; i<send_tmp.size(); i++)
        send_tmp[i] = (int) dmd.entrecv[i];
    getdatarecv(v, send_tmp, dmd.nbsd);
    find_datasend(dmd.entsend, v, dmd.entpart, dmd.nbsd, dmd.my_rank);
    for (Int n=0; n<nnbsd; n++) {
//         nbdmd[n].entrecv.resize(v[n].size());
//         for (Int i=0; i<v[n].size(); i++)
//             nbdmd[n].entrecv[i] = (Int) v[n][i];
        v[n].clear();
    }
    
    send_tmp.resize(dmd.vecrecv.size());
    for (Int i=0; i<send_tmp.size(); i++)
        send_tmp[i] = (int) dmd.vecrecv[i];
    getdatarecv(v, send_tmp, dmd.nbsd);
    find_datasend(dmd.vecsend, v, dmd.entpart, dmd.nbsd, dmd.my_rank);
    for (Int n=0; n<nnbsd; n++) {
//         nbdmd[n].vecrecv.resize(v[n].size());
//         for (Int i=0; i<v[n].size(); i++)
//             nbdmd[n].vecrecv[i] = (Int) v[n][i];
        v[n].clear();
    }
    
    send_tmp.resize(dmd.matrecv.size());
    for (Int i=0; i<send_tmp.size(); i++)
        send_tmp[i] = (int) dmd.matrecv[i];
    getdatarecv(v, send_tmp, dmd.nbsd);
    find_datasend(dmd.matsend, v, dmd.entpart, dmd.nbsd, dmd.my_rank);
    for (Int n=0; n<nnbsd; n++) {
//         nbdmd[n].matrecv.resize(v[n].size());
//         for (Int i=0; i<v[n].size(); i++)
//             nbdmd[n].matrecv[i] = (Int) v[n][i];
        v[n].clear();
    }
    
//     send_tmp.resize(dmd.ent2ind.size());
//     for (Int i=0; i<send_tmp.size(); i++)
//         send_tmp[i] = (int) dmd.ent2ind[i];
//     getdatarecv(v, send_tmp, dmd.nbsd);
//     for (Int n=0; n<nnbsd; n++) {
// //         nbdmd[n].ent2ind.resize(v[n].size());
// //         for (Int i=0; i<v[n].size(); i++)
// //             nbdmd[n].ent2ind[i] = (Int) v[n][i];
//         v[n].clear();
//     }
//     
//     send_tmp.resize(dmd.matsend.size());
//     for (Int i=0; i<send_tmp.size(); i++)
//         send_tmp[i] = (int) dmd.matsend[i];
//     getdatarecv(v, send_tmp, dmd.nbsd);
//     for (Int n=0; n<nnbsd; n++) {
// //         nbdmd[n].matsend.resize(v[n].size());
// //         for (Int i=0; i<v[n].size(); i++)
// //             nbdmd[n].matsend[i] = (Int) v[n][i];
//         v[n].clear();
//     }
//     
//     send_tmp.resize(dmd.entpart.size());
//     for (Int i=0; i<send_tmp.size(); i++)
//         send_tmp[i] = (int) dmd.entpart[i];
//     getdatarecv(v, send_tmp, dmd.nbsd);
//     for (Int n=0; n<nnbsd; n++) { 
//         nbdmd[n].entpart.resize(v[n].size());
//         for (Int i=0; i<v[n].size(); i++)
//             nbdmd[n].entpart[i] = (Int) v[n][i];
//         v[n].clear();
//     }
//     
//     send_tmp.resize(dmd.bcrs_rowent2ent.size());
//     for (Int i=0; i<send_tmp.size(); i++)
//         send_tmp[i] = (int) dmd.bcrs_rowent2ent[i];
//     getdatarecv(v, send_tmp, dmd.nbsd);
//     for (Int n=0; n<nnbsd; n++) {
//         nbdmd[n].bcrs_rowent2ent.resize(v[n].size());
//         for (Int i=0; i<v[n].size(); i++)
//             nbdmd[n].bcrs_rowent2ent[i] = (Int) v[n][i];
//         v[n].clear();
//     }
//     
//     send_tmp.resize(dmd.bcrs_colent2ent.size());
//     for (Int i=0; i<send_tmp.size(); i++)
//         send_tmp[i] = (int) dmd.bcrs_colent2ent[i];
//     getdatarecv(v, send_tmp, dmd.nbsd);
//     for (Int n=0; n<nnbsd; n++) {       
//         nbdmd[n].bcrs_colent2ent.resize(v[n].size());
//         for (Int i=0; i<v[n].size(); i++)
//             nbdmd[n].bcrs_colent2ent[i] = (Int) v[n][i];
//         v[n].clear();
//     }
}
        
void makedata(vector<Int> &datapts, vector<Int> &data, vector<Int> &nbsd)
{
    Int nnbsd = nbsd.size();   
    datapts.resize(nnbsd,0);
    
    Int i, n, k;    
    k = data.size();
    k = round(((double)k)/3);
    for (i=0; i<k; i++) 
        for (n = 0; n<nnbsd; n++)
            if (data[i] == nbsd[n])
                datapts[n] += 1;            
    
    // get only local numberings
    for (i=0; i<k; i++) 
        data[i] = data[k+i];
    data.resize(k);    
}

void sendrecvhdg(dmdstruct &dmd)
{      
    Int nent, j;
    
    nent = dmd.matrecv.size();    
    nent = round(((double) nent)/3);
    for (j=0; j<nent; j++)         
        dmd.matrecv[nent+j] = dmd.bcrs_rowent2ent[dmd.matrecv[nent+j]];           
        
    nent = dmd.matsend.size();    
    nent = round(((double) nent)/3);
    for (j=0; j<nent; j++)         
        dmd.matsend[nent+j] = dmd.bcrs_rowent2ent[dmd.matsend[nent+j]];     
    
    makedata(dmd.elemsendpts, dmd.elemsend, dmd.nbsd);
    makedata(dmd.elemrecvpts, dmd.elemrecv, dmd.nbsd);
    makedata(dmd.entsendpts, dmd.entsend, dmd.nbsd);
    makedata(dmd.entrecvpts, dmd.entrecv, dmd.nbsd);
    makedata(dmd.vecsendpts, dmd.vecsend, dmd.nbsd);
    makedata(dmd.vecrecvpts, dmd.vecrecv, dmd.nbsd);
    makedata(dmd.matsendpts, dmd.matsend, dmd.nbsd);
    makedata(dmd.matrecvpts, dmd.matrecv, dmd.nbsd);
}


void getmatrecv(vector< vector<int> > &nbmatrecv, vector<Int> &icpu, vector<Int> &jcpu, dmdstruct &dmd)
{        
    Int j, k, n, g, nent, nnbsd, ncpu;
    vector<Int> kcpu;
    
    nent = (Int) dmd.matrecv.size();    
    nent = round(((double) nent)/3);    
    nnbsd =  (Int) dmd.nbsd.size();
    icpu.resize(nent);
    jcpu.resize(nent);
    kcpu.resize(nnbsd);
    
    // temporary data
    vector< vector<int> > matrecv_tmp(nnbsd, vector<int>());    
    for (j=0; j<nnbsd; j++)
        kcpu[j] = 0;
    for (n=0; n<nent; n++) {
        ncpu = dmd.matrecv[0*nent+n]; // neighboring ncpu
        // find the index of ncpu
        for (k=0;k<nnbsd;k++)
            if (dmd.nbsd[k] == ncpu) {
                icpu[n] = k;       // index of ncpu        
                jcpu[n] = kcpu[k]; // count of ncpu
                break;
            } 
        kcpu[k] += 1;  
    }
    
    // allocate memory
    for (j=0; j<nnbsd; j++)
        matrecv_tmp[j].resize(kcpu[j]);
    
    for (n=0; n<nent; n++) {        
        g = dmd.matrecv[2*nent+n]; // global entity gj 
        k  = icpu[n];
        j  = jcpu[n];
        matrecv_tmp[k][j] = g;                                          
    }
    
//     if (dmd.my_rank==3) {        
//         cout<<dmd.my_rank<<"  "<<nent<<endl;
//         print2iarray(&dmd.matrecv[0], nent, 3);
//         for (j=0; j<nnbsd; j++)
//             print1iarray(&matrecv_tmp[j][0], (Int) matrecv_tmp[j].size());
//     }    

    getdatarecv(nbmatrecv, matrecv_tmp, dmd.nbsd);
}

void nbentcon(vector< vector<int> > &nbrowent2ent, vector< vector<int> > &nnrowent2ent, 
        vector< vector<int> > &nbentpart, vector< vector<int> > &nnentpart, 
        vector< vector<int> > &nbmatrecv, dmdstruct &dmd)
{    
    Int mj, nn, k, ng, gj, j, n, nk, jk, cjk, nnbsd;
    vector<Int> in;
    
    nnbsd =  (Int) dmd.nbsd.size();          
    vector< vector<int> > nbrowent2ent_tmp(nnbsd, vector<int>());    
    vector< vector<int> > nnrowent2ent_tmp(nnbsd, vector<int>());    
    vector< vector<int> > nbentpart_tmp(nnbsd, vector<int>());    
    vector< vector<int> > nnentpart_tmp(nnbsd, vector<int>());         
    
    // allocate memory
    for (n=0; n<nnbsd; n++) {        
        ng = nbmatrecv[n].size();    
        nbrowent2ent_tmp[n].resize(ng);
        nnrowent2ent_tmp[n].resize(ng);
        nnentpart_tmp[n].resize(ng);
        nk = 0;        
        for (j=0; j<ng; j++) {
            gj = nbmatrecv[n][j]; // global entity gj                              
            in = find(dmd.entpart,gj); // find index of global entity gj                              
            mj = in[0]; // local index of global entity gj
            nn = dmd.bcrs_rowent2ent[mj+1] - dmd.bcrs_rowent2ent[mj];    
            nk += nn;
        }        
        nbentpart_tmp[n].resize(nk);
    }    

    for (n=0; n<nnbsd; n++) { // for each neighboring cpu
        ng = nbmatrecv[n].size(); // number of global entities on neighboring cpu
        nk = 0;        
        for (j=0; j<ng; j++) { // for each global entity on neighboring cpu
            gj = nbmatrecv[n][j]; // global entity gj                              
            in = find(dmd.entpart,gj); // find index of global entity gj                              
            mj = in[0]; // local index of global entity gj
            nn = dmd.bcrs_rowent2ent[mj+1] - dmd.bcrs_rowent2ent[mj]; // number of entities connected to gj               
            nbrowent2ent_tmp[n][j] = dmd.bcrs_rowent2ent[mj]; // row index of gj
            nnrowent2ent_tmp[n][j] = nn; // number of entities connected to gj               
            nnentpart_tmp[n][j] = nk; //summing numbers of entities                
            for (k=0; k<nn; k++) {
                jk = dmd.bcrs_rowent2ent[mj]+k;
                // list of local neighboring entities to local entity mj   
                cjk = dmd.bcrs_colent2ent[jk]; 
                // list of global neighboring entities to gj on neighboring ncpu 
                //if (dmd.my_rank==0)
                //    cout<<n<<"  "<<nbentpart_tmp[n].size()<<"   "<<nk+k<<endl;
                nbentpart_tmp[n][nk+k] = dmd.entpart[cjk];                        
            }              
            nk += nn;
        }
    }    
    
    getdatarecv(nbrowent2ent, nbrowent2ent_tmp, dmd.nbsd);
    getdatarecv(nnrowent2ent, nnrowent2ent_tmp, dmd.nbsd);
    getdatarecv(nbentpart, nbentpart_tmp, dmd.nbsd);
    getdatarecv(nnentpart, nnentpart_tmp, dmd.nbsd);
}

void sendrecvedg(dmdstruct &dmd, vector< dmdstruct > &nbdmd)
{      
    Int j, k, m, nn, ij, nb, nj, mj, gj, gk, ik, ncpu, icpu, nent, irecv, nnbsd, nmax;            
    nnbsd =  (Int) dmd.nbsd.size();
    nent = (Int) dmd.matrecv.size();    
    nent = round(((double) nent)/3);
    
    k = (Int) dmd.bcrs_rowent2ent.size()-1;
    nmax = 0;
    for (j=0; j<k; j++) {
        nn = dmd.bcrs_rowent2ent[j+1] - dmd.bcrs_rowent2ent[j];
        nmax = (nmax>nn) ? nmax : nn;
    }    
    dmd.maxBlocksPerRow = nmax;
    
    //cout<<dmd.my_rank<<"  "<<nmax<<endl;
    //error("here");
            
    vector<Int> entnb, nbentj, nbelemj, nbelemk, nbelem, pj, gbent, mbent, in, rowj;
    vector<Int> isend(nnbsd,0);
    vector<Int> matrecv0(nent*nmax), matrecv1(nent*nmax), lmatrecv(nent);
    vector< vector<Int> > matsend0(nnbsd, vector<Int>());
    vector< vector<Int> > matsend1(nnbsd, vector<Int>());
    for (k=0;k<nnbsd;k++) {
        matsend0[k].resize(nent*nmax);
        matsend1[k].resize(nent*nmax);
    }
    for (j=0; j<nent; j++)
        lmatrecv[j] = dmd.matrecv[nent+j];
                
    vector< vector<int> > nbmatrecv(nnbsd, vector<int>());
    vector<Int> nbicpu;
    vector<Int> nbjcpu;
    getmatrecv(nbmatrecv, nbicpu, nbjcpu, dmd);
    vector< vector<int> > nbrowent2ent(nnbsd, vector<int>());
    vector< vector<int> > nnrowent2ent(nnbsd, vector<int>());
    vector< vector<int> > nbentpart(nnbsd, vector<int>());
    vector< vector<int> > nnentpart(nnbsd, vector<int>());        
    nbentcon(nbrowent2ent, nnrowent2ent, nbentpart, nnentpart, nbmatrecv, dmd);
    
    irecv = 0;
    for (j=0; j<nent; j++) {
        ncpu = dmd.matrecv[0*nent+j]; // neighboring cpu
        nj = dmd.matrecv[1*nent+j]; // local entity nj          
        gj = dmd.matrecv[2*nent+j]; // global entity gj                                  
        nn = dmd.bcrs_rowent2ent[nj+1] - dmd.bcrs_rowent2ent[nj];
        entnb.resize(nn);
        for (k=0; k<nn; k++) // for each entity neighboring the entity nj        
            entnb[k] = dmd.bcrs_colent2ent[dmd.bcrs_rowent2ent[nj]+k]; // entities neighboring the entity nj                        
        // take neighbors in lmatrecv
        nbentj = setintersection(entnb, lmatrecv);         
        // reorder the list of neighboring entities
        nn = (Int) nbentj.size();
        for (k=0; k<nn; k++)
            if (nbentj[k] == nj) {                
                nbentj.erase(nbentj.begin()+k);
                nbentj.insert(nbentj.begin(),nj);
                break;
            }                
        // get a subset of rj
        rowj.resize(nn);
        for (k=0; k<nn; k++) {
            in = find(entnb,nbentj[k]);            
            rowj[k] = dmd.bcrs_rowent2ent[nj]+in[0];            
        }        
        
//         entnb.insert(entnb.begin(),dmd.my_rank);
//         printiarray(entnb);
//         nbentj.insert(nbentj.begin(),dmd.my_rank);
//         printiarray(nbentj);
//         rowj.insert(rowj.begin(),dmd.my_rank);
//         printiarray(rowj);
//         error("here");
        
        // add [ncpu rowj(1)-1] to matrecv
        matrecv0[irecv] = ncpu;
        matrecv1[irecv] = rowj[0];
        irecv += 1;                                
       
        // find neighboring elements of the entity
        ij = dmd.ent2ind[gj];
        nn = dmd.rowent2elem[ij+1] - dmd.rowent2elem[ij];
        nbelemj.resize(nn);        
        // list of global neighboring elements to gj
        for (k=0; k<nn; k++) // for each elemennt neighboring the entity gj        
            nbelemj[k] = dmd.colent2elem[dmd.rowent2elem[ij]+k]; // elements neighboring the entity gj                
                        
        //cout<<dmd.my_rank<<"   "<<dmd.bcrs_rowent2elem[ij]<<"   "<<dmd.bcrs_rowent2elem[ij+1]<<endl;
        //cout<<dmd.my_rank<<"   "<<ncpu<<"   "<<gj<<"   "<<ij<<"   "<<nn<<endl;        
        //nbelemj.insert(nbelemj.begin(),dmd.my_rank);
        //printiarray(nbelemj); 
        //dmd.rowent2elem.insert(dmd.rowent2elem.begin(),dmd.my_rank);
        //printiarray(dmd.rowent2elem); 
        //error("here");
        
        // find the index of ncpu
        for (k=0;k<nnbsd;k++)
            if (dmd.nbsd[k] == ncpu)
                icpu = k;
        
        //nbelemj.insert(nbelemj.begin(),dmd.my_rank);
        //printiarray(nbelemj);
        //cout<<dmd.my_rank<<"   "<<icpu<<endl;
        //error("here");
        
        // CORRECTED UP TO HERE
        
                
        // local entity mj on neighboring ncpu              
//         nn = round(((double) nbdmd[icpu].matsend.size())/3);
//         for (k=0; k<nn; k++)
//             if (nbdmd[icpu].matsend[2*nn+k] == gj) {                
//                 mj = nbdmd[icpu].matsend[nn+k];     
//                 break;
//             }
//         in = find(nbdmd[icpu].entpart,gj);                              
//         mj = in[0];
//         nn = nbdmd[icpu].bcrs_rowent2ent[mj+1] - nbdmd[icpu].bcrs_rowent2ent[mj];    
//         pj.resize(nn);
//         mbent.resize(nn);
//         gbent.resize(nn);
//         for (k=0; k<nn; k++) {
//             pj[k] = nbdmd[icpu].bcrs_rowent2ent[mj]+k;
//             // list of local neighboring entities to local entity mj   
//             mbent[k] = nbdmd[icpu].bcrs_colent2ent[pj[k]]; 
//             // list of global neighboring entities to gj on neighboring ncpu 
//             gbent[k] = nbdmd[icpu].entpart[mbent[k]];      
//             
//             //cout<<gbent[k]<<"  "<<nbentpart[nbicpu[j]][nnentpart[nbicpu[j]][nbjcpu[j]]+k]<<endl;
//         }                                      
//         //cout<<nn<<"  "<<nnrowent2ent[nbicpu[j]][nbjcpu[j]]<<"  "<<nbdmd[icpu].bcrs_rowent2ent[mj]<<"  "<<nbrowent2ent[nbicpu[j]][nbjcpu[j]]<<endl;
        
        nn = nnrowent2ent[nbicpu[j]][nbjcpu[j]];
        nb = nbrowent2ent[nbicpu[j]][nbjcpu[j]];
        nj = nnentpart[nbicpu[j]][nbjcpu[j]];
        pj.resize(nn);        
        gbent.resize(nn);
        for (k=0; k<nn; k++) {
            pj[k] = nb+k;
            gbent[k] = nbentpart[nbicpu[j]][nj+k];            
        }
                
        //mbent = dmd{ncpu}.cbsr_colind(pj); 
        //cout<<dmd.my_rank<<"   "<<icpu<<"   "<<gj<<"   "<<in[0]<<"   "<<nn<<endl;
        //printiarray(nbdmd[icpu].bcrs_colent2ent);
                
        //pj.insert(pj.begin(),dmd.my_rank);
        //mbent.insert(mbent.begin(),dmd.my_rank);
        //gbent.insert(gbent.begin(),dmd.my_rank);
        //printiarray(pj);
        //printiarray(mbent);
        //printiarray(gbent);        
        //error("here");
        
        // add [i pj[0]] to matsend on neighboring ncpu                      
        matsend0[icpu][isend[icpu]] = dmd.my_rank;
        matsend1[icpu][isend[icpu]] = pj[0];
        isend[icpu] += 1;                        
        
        nb = (Int) rowj.size();
        for (k=1; k<nb; k++) {
            // global entity gk on self-cpu
            gk = dmd.entpart[nbentj[k]];               
            // local index of gk
            ik = dmd.ent2ind[gk]; 
            nn = dmd.rowent2elem[ik+1] - dmd.rowent2elem[ik];
            
//             if (dmd.my_rank==0)
//                 cout<<dmd.my_rank<<"   "<<gk<<"   "<<ik<<"   "<<nn<<"  "<<dmd.rowent2elem[ik]<<"  "<<dmd.rowent2elem[ik+1]<<endl;
                    
            // list of global neighboring elements to gk
            nbelemk.resize(nn);
            for (m=0; m<nn; m++) // for each elemennt neighboring the entity gj        
                nbelemk[m] = dmd.colent2elem[dmd.rowent2elem[ik]+m]; // elements neighboring the entity gj                
                        
            // intersect nbelemj and nbelemk to get a common set
            nbelem = setintersection(nbelemj, nbelemk);                     
            
//             if (dmd.my_rank==0) {
//                 printiarray(nbelem);
//                 printiarray(dmd.elempart);
//             }
            
            if (IsSubset(nbelem, dmd.elempart, (Int) nbelem.size(), (Int) dmd.elempart.size())==0) {
                // add [ncpu rowj(k)] to matrecv
                matrecv0[irecv] = ncpu;
                matrecv1[irecv] = rowj[k];
                irecv += 1;                        
                
                // find an element in gbent to match gk                
                in = find(gbent,gk);                
                // add [i pj(nk)] to matsend
                matsend0[icpu][isend[icpu]] = dmd.my_rank;
                matsend1[icpu][isend[icpu]] = pj[in[0]];
                isend[icpu] += 1;                                
            }            
        }
    }   
        
    dmd.matrecv.resize(irecv*3);
    for (j=0; j<irecv; j++) {
        dmd.matrecv[j] = matrecv0[j];
        dmd.matrecv[irecv+j] = matrecv1[j];
    }    
    
    vector< vector<int> > matsend(nnbsd, vector<int>());
    vector< vector<int> > matsend_(nnbsd, vector<int>());
    for (k=0;k<nnbsd;k++) {
        nn = isend[k];
        matsend_[k].resize(nn*3);
        for (j=0; j<nn; j++) {
            matsend_[k][j] = (int) matsend0[k][j];
            matsend_[k][nn+j] = (int) matsend1[k][j];
        }
    }
    getdatarecv(matsend, matsend_, dmd.nbsd);
    matsend_.clear();
    
    nmax = 0;
    for (k=0;k<nnbsd;k++) 
        nmax = nmax + round(((double) matsend[k].size())/3);
    dmd.matsend.resize(nmax*3);
    m = 0;
    for (k=0;k<nnbsd;k++) {
        nn = round(((double) matsend[k].size())/3);
        for (j=0; j<nn; j++) {            
            dmd.matsend[m] = (Int) matsend[k][j];
            dmd.matsend[nmax+m] = (Int) matsend[k][nn+j];
            dmd.matsend[2*nmax+m] = (Int) matsend[k][2*nn+j];
            m = m+1;
        }                
    }
    
    makedata(dmd.elemsendpts, dmd.elemsend, dmd.nbsd);
    makedata(dmd.elemrecvpts, dmd.elemrecv, dmd.nbsd);
    makedata(dmd.entsendpts, dmd.entsend, dmd.nbsd);
    makedata(dmd.entrecvpts, dmd.entrecv, dmd.nbsd);
    makedata(dmd.vecsendpts, dmd.vecsend, dmd.nbsd);
    makedata(dmd.vecrecvpts, dmd.vecrecv, dmd.nbsd);    
    makedata(dmd.matrecvpts, dmd.matrecv, dmd.nbsd);
    makedata(dmd.matsendpts, dmd.matsend, dmd.nbsd);
}

#endif

#endif


