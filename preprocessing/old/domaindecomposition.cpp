#ifndef __DOMAINDECOMPOSITION
#define __DOMAINDECOMPOSITION

#include <vector>
#include <algorithm>
#include <numeric>

#include "domaindecomposition.h"
#include "mathutilities.cpp"
#include "ioutil.cpp"

void overlappingsubdomain(vector<Int> &extintelem,vector<Int> &intelem,vector<Int> &t2t,Int overlappinglevel, Int ne, Int nee)
{            
    Int nint = (Int) intelem.size(); 
    Int i, j, k, n;
    
    //cout<<ne << "   "<<nee<<"   "<<nint<<endl;
    
    // element-to-element connectivity in nonoverlapping subdomain 
    vector<Int> intt2t;
    intt2t.resize(nint*nee);    
    for (i=0; i<nint; i++)
        for (j=0; j<nee; j++)        
            intt2t[j*nint+i] = t2t[j*ne+intelem[i]];            
        
//     print1iarray(&intt2t[0], (Int) intt2t.size());
//     error("here");
    
    // list of all elements connected to elements in nonoverlapping subdomain Omega_i
    vector<Int> elem;
    
    // sort and remove duplications
    //unique_copy(intt2t.begin(),intt2t.end(),elem.begin());    
    uniqueiarray(elem,intt2t);   
    
    // remove -1 
    if (elem[0] == -1)
        elem.erase(elem.begin());
        
//     print1iarray(&elem[0], (Int) elem.size());
//     error("here");
    
    // exterior elements connected to elements in nonoverlapping subdomain
    vector<Int> extelem = setdifference(elem, intelem);
    
//     print1iarray(&extelem[0], (Int) extelem.size());
//     error("here");
     
    // add more exterior elements if overlappinglevel >= 1    
    for (k=0;k<overlappinglevel;k++) {        
        n = (Int) extelem.size();
        intt2t.resize(n*nee);    
        for (i=0; i<n; i++)
            for (j=0; j<nee; j++)        
                intt2t[j*n+i] = t2t[j*ne+extelem[i]];            
        
        // sort and remove duplications
        uniqueiarray(elem,intt2t);    
        // remove -1 
        if (elem[0] == -1)
            elem.erase(elem.begin());

        set_difference(elem.begin(), elem.end(), intelem.begin(), intelem.end(), 
                        inserter(extelem, extelem.end()));        
        //set_difference(elem.begin(), elem.end(), intelem.begin(), intelem.end(), 
        //                inserter(sdelem, sdelem.end()));        
        //extelem.insert(extelem.begin(), sdelem.begin(), sdelem.end());        
    }

    // sort and remove duplications 
    uniqueiarray(elem,extelem);        
    
    // exterior elements    
    extintelem = setdifference(elem, intelem);
    
    // add interior elements to extintelem
    extintelem.insert(extintelem.begin(), intelem.begin(), intelem.end());                
    
//    print1iarray(&extintelem[0], (Int) extintelem.size());    
}

//function [rowent2elem,colent2elem,rowent2ent,colent2ent,ent,ent2ind] = elcon2entconmpi(elem2ent,elem)

void elcon2entcon(vector<Int> &rowent2elem,vector<Int> &colent2elem,vector<Int> &rowent2ent,vector<Int> &colent2ent,
        vector<Int> &ent2ind, vector<Int> &elem2ent, vector<Int> &elem, Int np)
{
        
    vector<Int> ent; 
    // sort and remove duplications
    uniqueiarray(ent,elem2ent);    
    
    Int i, j, ind;
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
    elc.resize(np);
    for (i=0; i<ne; i++) { // for each element i
        for (j=0;j<np;j++) // entities on element i  
            elc[j] = elem2ent[i*np+j];
        uniqueiarray(k,elc); // remove duplicate entities on element i      
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
        for (j=0;j<np;j++) // entities on element i  
            elc[j] = elem2ent[i*np+j];
        uniqueiarray(k,elc); // remove duplicate entities on element i      
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
        elc.resize(ni*np);
        for (j=0; j<ni; j++) // for each element neighboring the entity ei
        {
            e = colent2elem[rowent2elem[i]+j]; // element neighboring the entity ei
            e = elem2ind[e];
            for (a=0;a<np;a++) // entities on element i  
                elc[j*np+a] = elem2ent[e*np+a];                        
        }
        uniqueiarray(k,elc); // remove duplicate entities
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
       vector<Int> &elem2ent, Int np, Int ne)
{
    vector<Int> ent; 
    // sort and remove duplications
    uniqueiarray(ent,elem2ent);    

    Int i, j, ind;
    Int ndof = (Int) ent.size();    
    Int entmax = (Int) ent.back();        
    
    // store number of neighboring elements for each entity    
    rowent2elem.resize(ndof,0);
    for (i=0; i<ndof; i++)
        rowent2elem[i] = 0;
    
    vector<Int> elc, k;
    elc.resize(np);
    for (i=0; i<ne; i++) { // for each element i
        for (j=0;j<np;j++) // entities on element i  
            elc[j] = elem2ent[i*np+j];
        uniqueiarray(k,elc); // remove duplicate entities on element i      
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
        for (j=0;j<np;j++) // entities on element i  
            elc[j] = elem2ent[i*np+j];
        uniqueiarray(k,elc); // remove duplicate entities on element i      
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
    colent2ent.resize(ndof*np*me,0);    
    for (i=0; i<ndof*np*me; i++)
        colent2ent[i] = 0;
    Int ei, l, a, e, m = 0;
    for (i=0; i<ndof; i++) {   
        ei = ent[i];  // entity ei
        ni = rowent2elem[i+1] - rowent2elem[i];
        elc.resize(ni*np);
        for (j=0; j<ni; j++) // for each element neighboring the entity ei
        {
            e = colent2elem[rowent2elem[i]+j]; // element neighboring the entity ei
            for (a=0;a<np;a++) // entities on element i  
                elc[j*np+a] = elem2ent[e*np+a];                        
        }
        uniqueiarray(k,elc); // remove duplicate entities
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

void entelempart(vector<Int> &intelem, vector<Int> &intent, vector<Int> &rowent2elem, vector<Int> &colent2elem,
      vector<Int> &rowent2ent,vector<Int> &colent2ent, vector<Int> &ent2ind, 
      vector<Int> & elem2ent, vector<Int> & t2t, vector<Int> & elemall, vector<Int> & entall, 
      Int np, Int ne, Int nee, Int overlappinglevel, Int my_rank)
{
    // list of  elements in nonoverlapping subdomain my_rank                     
    intelem = find(elemall,my_rank); 
    
    // list of  entities in nonoverlapping subdomain my_rank                       
    intent = find(entall,my_rank);
    
    // list of elements on overlapping subdomain my_rank
    vector<Int> extintelem;        
    overlappingsubdomain(extintelem,intelem,t2t,overlappinglevel,ne,nee);
    
    // list of entities on overlapping subdomain my_rank
    Int i, j, n1 = (Int) extintelem.size();
    vector<Int> extintelcon(n1*np,0); 
    for (i=0; i<n1; i++)
        for (j=0; j<np; j++)
            extintelcon[i*np+j] = elem2ent[extintelem[i]*np+j];    
    
//    print2iarray(&extintelcon[0], np, n1);    
//s    error("here");
    
    // (rowent2elem,colent2elem) global entity-to-element connectivity
    // (rowent2ent,colent2ent) global entity-to-entity connectivity  
    // ent2ind: entity-to-index mapping
    elcon2entcon(rowent2elem,colent2elem,rowent2ent,colent2ent,ent2ind,extintelcon,extintelem,np);    
}

void find_intfent(vector<Int> &intfent, vector<Int> &elemintfent, vector<Int> &intent, 
        vector<Int> &ent2ind, vector<Int> &rowent2elem, vector<Int> &colent2elem,
        vector<Int> &rowent2ent, vector<Int> &colent2ent)
{
    //find all interface entities in intent
    Int i, j, nj, ij, ni, m, nent;
    
    nent = intent.size();
    intfent.resize(nent);    
    for (i=0; i<nent; i++)
        intfent[i] = 0;    
    
    vector<Int> entnb;        
    m = 0; 
    for(i=0;i<nent;i++) { // For each entity in intent
        nj = intent[i]; // entity nj
        ij = ent2ind[nj]; // index of nj
        ni = rowent2ent[ij+1] - rowent2ent[ij];
        entnb.resize(ni); // list of neirghboring entities 
        for (j=0; j<ni; j++) // for each entity neighboring the entity ei        
            entnb[j] = colent2ent[rowent2ent[ij]+j]; // entity neighboring the entity ei                                
        if (IsSubset(entnb, intent) == 0) { // entnb is not a subset of intent
            intfent[m] = intent[i];  
            m += 1;
        }        
    }
    intfent.resize(m);    

   // cout<<nent<<"  "<<m<<endl;    
    //print1iarray(&intent[0], (Int) intent.size());
    //print1iarray(&intfent[0], (Int) intfent.size());
    //print1iarray(&elemintfent[0], (Int) elemintfent.size());
        
    vector<Int> elem;
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
    }    
    // sort and remove duplications 
    uniqueiarray(elemintfent,elem); // elements connected to interface entities           
}

void find_extelem(vector<Int> &intfelem, vector<Int> &extelem, vector<Int> &elemintfent, 
        vector<Int> &intelem, vector<Int> &t2t, Int preconditioner, Int ne, Int nee)
{
    //interface elements
    intfelem = setintersection(elemintfent, intelem);                     
    
    // exterior elements
    extelem = setdifference(elemintfent, intfelem);    
        
    Int ni = intfelem.size();
    if (preconditioner>0) {
        Int i, j, k;
        
        vector<Int> el(ni*nee);
        for (i=0; i<ni; i++)
            for (j=0; j<nee; j++)
                el[i*nee+j] = t2t[j*ne+intfelem[i]];
        
        // sort and remove duplications 
        vector<Int> elem;
        uniqueiarray(elem,el);        
        if (elem[0]<0)
            elem.erase(elem.begin());
        
        // exterior elements
        set_difference(intelem.begin(), intelem.end(), 
                       elem.begin(), elem.end(), 
                       inserter(extelem, extelem.end()));    
        
        for (k=1; k<preconditioner; k++) {
            ni = extelem.size();
            el.resize(ni*nee);
            for (i=0; i<ni; i++)
                for (j=0; j<nee; j++)
                    el[i*nee+j] = t2t[j*ne+intfelem[i]];
            
            // sort and remove duplications 
            uniqueiarray(elem,el);        
            if (elem[0]<0)
                elem.erase(elem.begin());

            // exterior elements
            set_difference(intelem.begin(), intelem.end(), 
                           elem.begin(), elem.end(), 
                           inserter(extelem, extelem.end()));                
        }        
    }
    extelem = uniqueiarray(extelem);
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
                  vector<Int> &elemall)
{        
    // store elements received from neighboring subdomains to assemble linear system
    Int i, k, ei;
    Int nc = extelem.size();
    Int nb = intelem.size();
    elemrecv.resize(nc*3);
    for (i=0; i<nc; i++) {
        ei = extelem[i];  // global entity numbering
        k = elemall[ei]; // neighboring subdomain k
        elemrecv[0*nc+i] = k;  // store neighboring subdomain k
        elemrecv[1*nc+i] = nb+i;  // store local entity numbering
        elemrecv[2*nc+i] = ei;  // store global entity numbering
    }                
    //sort_columnk(elemrecv, nc, 3, 0);    
    sort_rows(elemrecv, nc, 3);
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

void find_entrecv(vector<Int> &entrecv, vector<Int> &nbsd, vector<Int> &extent, vector<Int> &intent,
                  vector<Int> &entall)
{    
    // store entities received from neighboring subdomains
    Int i, k, ei;
    Int nc = extent.size();
    Int nb = intent.size();
    entrecv.resize(nc*3);
    vector<Int> entrecv0(nc,0);
    for (i=0; i<nc; i++) {
        ei = extent[i];  // global entity numbering
        k = entall[ei]; // neighboring subdomain k
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
    vector<Int> elemnb, on(na,0); 
    for (i=0; i<na; i++) {
        nj = entrecv[2*na+i]; // exterior global entity nj
        ij = ent2ind[nj]; // index of nj (local entity)
        ni = rowent2elem[ij+1] - rowent2elem[ij];
        elemnb.resize(ni); // list of neirghboring elements 
        for (j=0; j<ni; j++) // for each element neighboring the entity ei        
            elemnb[j] = colent2elem[rowent2elem[ij]+j]; // entity neighboring the entity ei                         
        // entity nj is connected to interface entities
        if (IsSubset(elemnb,elempart)==0) {
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
         vector<Int> &entpart, vector<Int> &elempart, vector<Int> &elconhdg, vector<Int> &t2fhdg, Int nfe, Int npf)
{
    Int i, j, k, ent;
    Int ne = elempart.size();   
    Int nn = round(((double) t2fhdg.size()/nfe));
    
    // local t2f for hdg
    t2f.resize(ne*nfe);
    for(j=0;j<nfe;j++) 
        for (i=0;i<ne;i++) {
            ent = t2fhdg[j*nn+elempart[i]];
            for (k=0; k<entpart.size(); k++)
                if (entpart[k]==ent) {
                    t2f[j*ne+i] = k;  
                    break;
                }
        }        
         
    // local elcon for hdg
    vector<Int> ii(npf), ind(npf);
    elcon.resize(npf*nfe*ne);    
    for(i=0;i<ne;i++) 
        for (j=0;j<nfe;j++) {
            for (k=0;k<npf;k++) {
                ii[k] = elconhdg[elempart[i]*nfe*npf+j*npf+k] - t2fhdg[j*nn+elempart[i]]*npf;
                ind[k] = t2f[j*ne+i]*npf+k;
            }
            for (k=0;k<npf;k++)
                elcon[i*nfe*npf+j*npf+k] = ind[ii[k]];
        }   
    
    // take transpose
    vector<Int> t2ft(ne*nfe);
    for (i=0;i<ne;i++)
        for(j=0;j<nfe;j++)         
            t2ft[i*nfe+j] = t2f[j*ne+i];
    
    // local entity-to-element and entity-to-entity connectivities
    elcon2entcon(bcrs_rowent2elem,bcrs_colent2elem,bcrs_rowent2ent,bcrs_colent2ent, t2ft, nfe, ne);                                
}

void make_elconedg(vector<Int> &elcon, vector<Int> &t2f, vector<Int> &bcrs_rowent2elem,
         vector<Int> &bcrs_colent2elem, vector<Int> &bcrs_rowent2ent, vector<Int> &bcrs_colent2ent,
         vector<Int> &entpart, vector<Int> &elempart, vector<Int> &elconedg, vector<Int> &t2fedg, Int nfe, Int npf)
{
    
    Int i, j;
    Int ne = (Int) elempart.size();      
    Int nn = round(((double) t2fedg.size()/nfe));
    t2f.resize(ne*nfe);
    for(j=0;j<nfe;j++) 
        for (i=0;i<ne;i++) 
            t2f[j*ne+i] = t2fedg[j*nn+elempart[i]];

    vector<Int> facesInProcessor;
    uniqueiarray(facesInProcessor,t2f); 
    
//      print2iarray(&t2f[0],ne,nfe);
//      printiarray(facesInProcessor);
    
    Int maxface = facesInProcessor.back();
    Int nface = (Int) facesInProcessor.size();
    vector<Int> facemapping(maxface+1,0);
    for (i=0;i<nface;i++)
        facemapping[facesInProcessor[i]] = i;
    for(j=0;j<nfe;j++) 
        for (i=0;i<ne;i++) 
            t2f[j*ne+i] = facemapping[t2f[j*ne+i]];
    
    //print2iarray(&t2f[0],ne,nfe);
        
    elcon.resize(npf*nfe*ne);    
    for(i=0;i<ne;i++) 
        for (j=0;j<nfe*npf;j++)             
            elcon[i*nfe*npf+j] = elconedg[elempart[i]*nfe*npf+j];      
    
    Int maxent = *max_element(entpart.begin(), entpart.end()); 
    Int nent = (Int) entpart.size();
    vector<Int> entmapping(maxent+1,0);
    for (i=0;i<nent;i++)
        entmapping[entpart[i]] = i;
    
    for(i=0;i<ne;i++) 
        for (j=0;j<nfe*npf;j++)             
            elcon[i*nfe*npf+j] = entmapping[elcon[i*nfe*npf+j]];
        
    //print2iarray(&elcon[0],npf*nfe,ne);
    
    // local entity-to-element and entity-to-entity connectivities    
    elcon2entcon(bcrs_rowent2elem,bcrs_colent2elem,bcrs_rowent2ent,bcrs_colent2ent, elcon, npf*nfe, ne);                                
}

void find_datasend(vector<Int> &datasend,vector< vector<Int> > &datarecv,vector<Int> &epart,vector<Int> &nbsd,Int my_rank)
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
        nrecv = datarecv[n].size();
        nrows = round(((double)nrecv)/3);         
        for (i=0; i<nrows; i++) {            
            if (datarecv[n][i]==my_rank) {
                e = datarecv[n][2*nrows+i];
                k = find(epart,e)[0];
                datasend[0*ne+j] = nbsd[n];
                datasend[1*ne+j] = k;
                datasend[2*ne+j] = epart[k];
                j = j + 1;
            }
        }            
    }              
}

void domaindecomposition(dmdstruct &dmd, vector<Int> &elem2ent, vector<Int> &elconhdg, 
        vector<Int> &t2t, vector<Int> &t2f, vector<Int> &elemall, vector<Int> &entall, 
        Int preconditioner, Int overlappinglevel, Int hdg, Int my_rank)
{
   
    Int ne = elemall.size();
    Int nt2t = t2t.size();
    Int nee = round(((double) nt2t)/ne);
    Int nel = elem2ent.size();
    Int np = round(((double)nel)/ne);      
    Int nt2f = t2f.size();
    Int nfe = nt2f/ne;
    Int net = elconhdg.size();
    Int npf = round(((double) net)/(ne*nfe));
    
    vector<Int> intfent, intfelem, extelem, extent, elemintfent;
   
    dmd.my_rank = my_rank;
    
    entelempart(dmd.intelem, dmd.intent, dmd.rowent2elem, dmd.colent2elem, dmd.rowent2ent, dmd.colent2ent, 
               dmd.ent2ind, elem2ent, t2t, elemall, entall, np, ne, nee, overlappinglevel, my_rank);

    find_intfent(intfent, elemintfent, dmd.intent, dmd.ent2ind, dmd.rowent2elem, dmd.colent2elem,
                dmd.rowent2ent, dmd.colent2ent);

    find_extelem(intfelem, extelem, elemintfent, dmd.intelem, t2t, preconditioner, ne, nee);
  
    make_elempart(dmd.elempart, dmd.elempartpts, dmd.intelem, intfelem, extelem);      
   
    find_extent(extent, intfelem, extelem, dmd.intent, elem2ent, np);
   
    make_elempart(dmd.entpart, dmd.entpartpts, dmd.intent, intfent, extent);
   
    find_elemrecv(dmd.elemrecv, extelem, dmd.intelem, elemall);
   
    find_entrecv(dmd.entrecv, dmd.nbsd, extent, dmd.intent, entall);
   
    find_vecrecv(dmd.vecrecv, extent, intfent, dmd.entrecv, dmd.ent2ind, dmd.rowent2ent, dmd.colent2ent);
   
    find_matrecv(dmd.matrecv, extent, dmd.elempart, dmd.entrecv, dmd.ent2ind, dmd.rowent2elem, dmd.colent2elem);   
    
    if (hdg==1)
        make_elconhdg(dmd.elcon, dmd.t2f, dmd.bcrs_rowent2elem, dmd.bcrs_colent2elem, 
            dmd.bcrs_rowent2ent, dmd.bcrs_colent2ent, dmd.entpart, dmd.elempart, 
            elconhdg, t2f, nfe, npf);
    else
        make_elconedg(dmd.elcon, dmd.t2f, dmd.bcrs_rowent2elem, dmd.bcrs_colent2elem, 
            dmd.bcrs_rowent2ent, dmd.bcrs_colent2ent, dmd.entpart, dmd.elempart, 
            elem2ent, t2f, nfe, npf);
};

void domaindecompositionnb(dmdstruct &dmd, vector< dmdstruct > &nbdmd, vector<Int> &elem2ent, 
        vector<Int> &elconhdg,  vector<Int> &t2t, vector<Int> &t2f, vector<Int> &elemall, 
        vector<Int> &entall, Int preconditioner, Int overlappinglevel, Int hdg)
{
   
    // TO DO: Use MPI instead
    Int i, nnbsd = (Int) dmd.nbsd.size();
    for (i=0; i< nnbsd; i++) // domain decomposition for neighboring cpus       
        domaindecomposition(nbdmd[i], elem2ent, elconhdg, t2t, t2f, elemall, entall, 
                preconditioner, overlappinglevel, hdg, dmd.nbsd[i]);      
    
    vector< vector<Int> > nbrecv;
    
    // data to send to neighboring cpus       
    for (i=0; i< nnbsd; i++)     
        nbrecv[i] = nbdmd[i].elemrecv;    
    find_datasend(dmd.elemsend, nbrecv, dmd.elempart, dmd.nbsd, dmd.my_rank);
    
    for (i=0; i< nnbsd; i++)     
        nbrecv[i] = nbdmd[i].entrecv;    
    find_datasend(dmd.entsend, nbrecv, dmd.entpart, dmd.nbsd, dmd.my_rank);
    
    for (i=0; i< nnbsd; i++)     
        nbrecv[i] = nbdmd[i].vecrecv;    
    find_datasend(dmd.vecsend, nbrecv, dmd.entpart, dmd.nbsd, dmd.my_rank);
    
    for (i=0; i< nnbsd; i++)     
        nbrecv[i] = nbdmd[i].matrecv;    
    find_datasend(dmd.matsend, nbrecv, dmd.entpart, dmd.nbsd, dmd.my_rank);    
};

#ifdef HAVE_MPI
vector< vector<Int> > getdatarecv(vector<Int> & datasend, vector<Int> & nbsd)
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
    
    /* Non-blocking send */    
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];   
        if (nsend>0) {
            MPI_Isend(&datasend[0], nsend, MPI_INT, neighbor, 0,
                   MPI_COMM_WORLD, &requests[requestCounter]);
            requestCounter += 1;
        }        
    }            
    
    // initialize datarecv for all neighboring cpus   
    vector< vector<Int> > datarecv(nnbsd, vector<Int>());            
    
    /* Non-blocking receive */
    Int nrecv;
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];
        
        // Probe for an incoming message from neighboring cpu
        MPI_Probe(neighbor, 0, MPI_COMM_WORLD, &statuses[requestCounter]);
        
        // When probe returns, the status object has the size and other
        // attributes of the incoming message. Get the message size
        MPI_Get_count(&statuses[requestCounter], MPI_INT, &nrecv);                

        if (nrecv>0) {
            datarecv[n].resize(nrecv);
            MPI_Irecv(&datarecv[n][0], nrecv, MPI_INT, neighbor, 0,
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

vector< vector<Int> > getdatarecv(vector < vector<Int> > & datasend, vector<Int> & nbsd)
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
    
    /* Non-blocking send */    
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];   
        nsend = (Int) datasend[n].size();
        if (nsend>0) {
            MPI_Isend(&datasend[n][0], nsend, MPI_INT, neighbor, 0,
                   MPI_COMM_WORLD, &requests[requestCounter]);
            requestCounter += 1;
        }        
    }            
    
    // initialize datarecv for all neighboring cpus   
    vector< vector<Int> > datarecv(nnbsd, vector<Int>());            
    
    /* Non-blocking receive */
    Int nrecv;
    for (n=0; n<nnbsd; n++) {
        neighbor = nbsd[n];
        
        // Probe for an incoming message from neighboring cpu
        MPI_Probe(neighbor, 0, MPI_COMM_WORLD, &statuses[requestCounter]);
        
        // When probe returns, the status object has the size and other
        // attributes of the incoming message. Get the message size
        MPI_Get_count(&statuses[requestCounter], MPI_INT, &nrecv);                

        if (nrecv>0) {
            datarecv[n].resize(nrecv);
            MPI_Irecv(&datarecv[n][0], nrecv, MPI_INT, neighbor, 0,
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

void make_nbdmd(vector< dmdstruct > &nbdmd, dmdstruct &dmd)
{
    Int nnbsd = (Int) dmd.nbsd.size();                        
    vector< vector<Int> > v;
    
    v = getdatarecv(dmd.elemrecv, dmd.nbsd);            
    find_datasend(dmd.elemsend, v, dmd.elempart, dmd.nbsd, dmd.my_rank);
    for (Int n=0; n<nnbsd; n++) 
        nbdmd[n].elemrecv = v[n];        
    
    v = getdatarecv(dmd.entrecv, dmd.nbsd);
    find_datasend(dmd.entsend, v, dmd.entpart, dmd.nbsd, dmd.my_rank);
    for (Int n=0; n<nnbsd; n++) 
        nbdmd[n].entrecv = v[n];
    
    v = getdatarecv(dmd.vecrecv, dmd.nbsd);
    find_datasend(dmd.vecsend, v, dmd.entpart, dmd.nbsd, dmd.my_rank);
    for (Int n=0; n<nnbsd; n++) 
        nbdmd[n].vecrecv = v[n];
    
    v = getdatarecv(dmd.matrecv, dmd.nbsd);
    find_datasend(dmd.matsend, v, dmd.entpart, dmd.nbsd, dmd.my_rank);
    for (Int n=0; n<nnbsd; n++) 
        nbdmd[n].matrecv = v[n];
                
    v = getdatarecv(dmd.ent2ind, dmd.nbsd);
    for (Int n=0; n<nnbsd; n++) 
        nbdmd[n].ent2ind = v[n];                
    
    v = getdatarecv(dmd.matsend, dmd.nbsd);
    for (Int n=0; n<nnbsd; n++) 
        nbdmd[n].matsend = v[n];
    
    v = getdatarecv(dmd.entpart, dmd.nbsd);
    for (Int n=0; n<nnbsd; n++) 
        nbdmd[n].entpart = v[n];                
    
    v = getdatarecv(dmd.bcrs_rowent2ent, dmd.nbsd);
    for (Int n=0; n<nnbsd; n++) 
        nbdmd[n].bcrs_rowent2ent = v[n];
    
    v = getdatarecv(dmd.bcrs_colent2ent, dmd.nbsd);
    for (Int n=0; n<nnbsd; n++)         
        nbdmd[n].bcrs_colent2ent = v[n];              
}

#endif
        
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
        in = find(nbdmd[icpu].entpart,gj);                              
        mj = in[0];
        nn = nbdmd[icpu].bcrs_rowent2ent[mj+1] - nbdmd[icpu].bcrs_rowent2ent[mj];    
        pj.resize(nn);
        mbent.resize(nn);
        gbent.resize(nn);
        for (k=0; k<nn; k++) {
            pj[k] = nbdmd[icpu].bcrs_rowent2ent[mj]+k;
            // list of local neighboring entities to local entity mj   
            mbent[k] = nbdmd[icpu].bcrs_colent2ent[pj[k]]; 
            // list of global neighboring entities to gj on neighboring ncpu 
            gbent[k] = nbdmd[icpu].entpart[mbent[k]];
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
            
            if (IsSubset(nbelem, dmd.elempart)==0) {
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
        
    vector< vector<Int> > matsend(nnbsd, vector<Int>());
    for (k=0;k<nnbsd;k++) {
        nn = isend[k];
        matsend[k].resize(nn*3);
        for (j=0; j<nn; j++) {
            matsend[k][j] = matsend0[k][j];
            matsend[k][nn+j] = matsend1[k][j];
        }
    }        
    matsend = getdatarecv(matsend, dmd.nbsd);
    
    nmax = 0;
    for (k=0;k<nnbsd;k++) 
        nmax = nmax + round(((double) matsend[k].size())/3);
    dmd.matsend.resize(nmax*3);
    m = 0;
    for (k=0;k<nnbsd;k++) {
        nn = round(((double) matsend[k].size())/3);
        for (j=0; j<nn; j++) {            
            dmd.matsend[m] = matsend[k][j];
            dmd.matsend[nmax+m] = matsend[k][nn+j];
            dmd.matsend[2*nmax+m] = matsend[k][2*nn+j];
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

// void (dmdstruct &dmd, vector< dmdstruct > &nbdmd, Int my_rank)
// {
//     find_datasend(dmd.elemsend, nbelemrecv, dmd.elempart, dmd.nbsd, my_rank);
// }


// void entelempart(vector<Int> &intelem, vector<Int> &intent, vector<Int> &rowent2elem, vector<Int> &colent2elem,
//       vector<Int> &rowent2ent,vector<Int> &colent2ent, vector<Int> &ent2ind, 
//       vector<Int> & elcon, vector<Int> & t2t, vector<Int> & elemall, vector<Int> & entall, 
//       Int np, Int ne, Int nfe, Int overlappinglevel, Int my_rank)

      
// void domaindecompositionmpi(elcon,t2t,t2f,nproc,preconditioner,elconhdg,my_rank)
// {        
//     
//     // list of  elements in nonoverlapping subdomain my_rank                     
//     Int nelem = elemall.size();    
//     vector<Int> intelem(nelem,0);
//     j = 0;
//     for (i=0; i<nelem; i++) 
//         if (elemall[i] = my_rank) {
//             intelem[j] = i;
//             j += 1;
//         }    
//     intelem.erase(intelem.begin()+j,intelem.end());
//     
//     // list of  entities in nonoverlapping subdomain my_rank                       
//     Int nent = entall.size();
//     vector<Int> intent(nent,0);
//     j = 0;
//     for (i=0; i<nent; i++) 
//         if (entall[i] = my_rank) {
//             intent[j] = i;
//             j += 1;
//         }
//     intent.erase(intent.begin()+j,intent.end());
//     
//     vector<Int> tmp;
//     Int nintelem = intelem.size();
//     for (i=0; i<nintelem; i++) { 
//         m = intelem[i]*np;
//         for (j=0; j<np; j++)
//             elc[j] = elcon[m+j];
//         sort(elc.begin(),elc.end());   
//         set_intersection(elc.begin(), elc.end(), intent.begin(), intent.end(), 
//                         inserter(tmp, tmp.begin()));    
//          if (tmp.empty()) {
//             cout<<"Element: "<<i<<"\n"; 
//             error("An element in the processor did not have any entity in the processors. The issue need to be fixed.");        
//          }
//     }
//                 
//     // list of elements on overlapping subdomain i
//     vector<Int> extintelem;        
//     overlappingsubdomain(extintelem,intelem,t2t,overlappinglevel,ne,nfe)
//            
//     Int n1 = extintelem.size();
//     vector<Int> exelcon(n1*np,0); 
//     for (i=0; i<n1; i++)
//         for (j=0; j<np; j++)
//             exelcon[i*np+j] = elcon[extintelem[i]*np+j];    
//     // (rowent2elem,colent2elem) entity-to-element connectivity
//     // (rowent2ent,colent2ent) entity-to-entity connectivity        
//     elcon2entcon(rowent2elem,colent2elem,rowent2ent,colent2ent, ent2ind, exelcon, extintelem);
// 
//     
//     //find all interface entities in dmd{i}.intent
//     vector<Int> tmp, intfent(nent,0);        
//     vector<Int>::iterator it;
//     m = 0;
//     for(i=0;i<nent;i++) { // For each entity in intent
//         nj = intent[i]; // entity nj
//         ij = ent2ind[nj]; // index of nj
//         ni = rowent2elem[ij+1] - rowent2elem[ij];
//         entnb.resize(ni,0); // list of neirghboring entities 
//         for (j=0; j<ni; j++) // for each entity neighboring the entity ei        
//             entnb[j] = colent2elem[rowent2elem[ij]+j]; // entity neighboring the entity ei                
//         // get elements of entnb which are not in intent    
//         it = set_difference(entnb.begin(), entnb.end(), intent.begin(), intent.end(), tmp.begin());    
//         tmp.resize(it-tmp.begin());
//         
//         if (tmp.size() == 0) { // entity nj is on the interface between two subdomains                  
//             intfent[m] = intent[i];  
//             m += 1;
//         }
//     }
//     intfent.resize(m);    
//     
//     elem.resize(0);
//     for (i=0; i<m; i++) {
//         nj = intfent[i]; // entity nj
//         ij = ent2ind[nj]; // index of nj
//         ni = rowent2elem[ij+1] - rowent2elem[ij];        
//         for (j=0; j<ni; j++) // for each entity neighboring the entity ei        
//             entnb[j] = colent2elem[rowent2elem[ij]+j]; // entity neighboring the entity ei                
//         elem.insert(elem.end(), entnb.begin(), entnb.begin()+ni);
//     }    
//     // sort and remove duplications 
//     uniqueiarray(elemNeighboringInterfEnt,elem);        
//     
//     elem.resize(0);
//     for (i=0; i<m; i++) {
//         nj = intfent[i]; // entity nj
//         ij = ent2ind[nj]; // index of nj
//         ni = rowent2ent[ij+1] - rowent2ent[ij];
//         for (j=0; j<ni; j++) // for each entity neighboring the entity ei        
//             entnb[j] = colent2ent[rowent2ent[ij]+j]; // entity neighboring the entity ei                
//         elem.insert(elem.end(), entnb.begin(), entnb.begin()+ni);
//     }    
//     // sort and remove duplications 
//     uniqueiarray(entNeighboringInterfEnt,elem);        
//         
//     //interface elements
//     set_intersection(elemNeighboringInterfEnt.begin(), elemNeighboringInterfEnt.end(), 
//                      intelem.begin(), intelem.end(), 
//                      inserter(intfelem, intfelem.begin()));    
//     
//     // exterior elements
//     set_difference(elemNeighboringInterfEnt.begin(), elemNeighboringInterfEnt.end(), 
//                      intelem.begin(), intelem.end(), 
//                      inserter(extelem, extelem.begin()));    
//     
//     ni = intfelem.size();
//     if (preconditioner>0) {
//         vector<Int> el(ni*nb);
//         for (i=0; i<ni; i++)
//             for (j=0; j<nb; j++)
//                 el[i*nb+j] = t2t[j*ne+intfelem[i]];
//         // sort and remove duplications 
//         uniqueiarray(elem,el);        
//         if (elem[0]==0)
//             elem.erase(elem.begin());
//         
//         // exterior elements
//         set_difference(intelem.begin(), intelem.end(), 
//                        elem.begin(), elem.end(), 
//                        inserter(extelem, extelem.end()));    
//         
//         for (k=1; k<preconditioner; k++) {
//             ni = extelem.size();
//             el.resize(ni*nb);
//             for (i=0; i<ni; i++)
//                 for (j=0; j<nb; j++)
//                     el[i*nb+j] = t2t[j*ne+intfelem[i]];
//             // sort and remove duplications 
//             uniqueiarray(elem,el);        
//             if (elem[0]==0)
//                 elem.erase(elem.begin());
// 
//             // exterior elements
//             set_difference(intelem.begin(), intelem.end(), 
//                            elem.begin(), elem.end(), 
//                            inserter(extelem, extelem.end()));                
//         }        
//     }
//     unique(extelem.begin(),extelem.end());
// 
//     // [interior elements, interface elements, exterior elements]  
//     it = set_difference(intelem.begin(), intelem.end(), 
//                        intfelem.begin(), intfelem.end(), 
//                        inserter(elempart, elempart.begin()));                    
//     elempart.resize(it-elempart.begin());
//     elempartpts[0] = elempart.size();    
//     elempart.insert(elempart.end(), intfelem.begin(), intfelem.end());
//     elempartpts[1] = intfelem.size();
//     elempart.insert(elempart.end(), extelem.begin(), extelem.end());
//     elempartpts[2] = extelem.size();
//     
//     // store elements received from neighboring subdomains to assemble linear system
//     na = extelem.size();
//     nb = intelem.size();
//     vector<Int> elemrecv(na*3);
//     for (i=0; i<na; i++) {
//         ei = extelem[i];  // global element numbering
//         k = elemall[ei]; // neighboring subdomain k
//         elemrecv[0*na+i] = k;  // store neighboring subdomain k
//         elemrecv[1*na+i] = nb+i;  // store local element numbering
//         elemrecv[2*na+i] = ei;  // store global element numbering
//     }        
//    sort_columnk(elemrecv, na, 3, 0);
//         
//     // find all entities that are connected to interface elements and exterior elements    
//     nc = elempartpts[0];
//     nd = na+nb+nc;
//     vector<Int> ent((na+nb)*np);
//     for (i=nc;i<nd;i++) {
//         m = elempart[i]*np;        
//         for (j=0; j<np; j++)
//             ent[(i-nc)*np+j] = elcon[m+j];            
//     }
//     unique(ent.begin(),ent.end());
//     // exterior entities
//     it = set_difference(ent.begin(), ent.end(), 
//                         intent.begin(), intent.end(), 
//                         inserter(extent, extent.begin()));                    
//     extent.resize(it-extent.begin());
//     
//     // [interior entities, interface entities, exterior entities]
//     it = set_difference(intent.begin(), intent.end(), 
//                        intfent.begin(), intfent.end(), 
//                        inserter(entpart, entpart.begin()));                    
//     entpart.resize(it-entpart.begin());
//     entpartpts[0] = entpart.size();    
//     entpart.insert(entpart.end(), intfent.begin(), intfent.end());
//     entpartpts[1] = intfent.size();
//     entpart.insert(entpart.end(), extent.begin(), extent.end());
//     entpartpts[2] = extent.size();
//     
//     // store entities received from neighboring subdomains
//     na = extent.size();
//     nb = intent.size();
//     vector<Int> entrecv(na*3);
//     for (i=0; i<na; i++) {
//         ei = extent[i];  // global entity numbering
//         k = entall[ei]; // neighboring subdomain k
//         entrecv[0*na+i] = k;  // store neighboring subdomain k
//         entrecv[1*na+i] = nb+i;  // store local entity numbering
//         entrecv[2*na+i] = ei;  // store global entity numbering
//     }                
//     sort_columnk(entrecv, na, 3, 0);
//     // list of neighboring subdomains from which info is required for matrix-vector product    
//     unique_copy(entrecv.begin(),entrecv.begin()+na,nbsd); 
//     
//     // find exterior entities connected to interface entities to perform matrix-vector product. 
//     vector<Int> on(na,0); 
//     k = 0;
//     for (i=0; i<na; i++) {
//         nj = extent[i]; // exterior entity nj
//         ij = ent2ind[nj]; // index of nj
//         ni = rowent2ent[ij+1] - rowent2ent[ij];
//         for (j=0; j<ni; j++) // for each entity neighboring the entity ei        
//             entnb[j] = colent2ent[rowent2ent[ij]+j]; // entity neighboring the entity ei                        
//         it = set_intersection(entnb.begin(), entnb.end(), intfent.begin(), intfent.end(), 
//                         inserter(ent, ent.begin()));    
//         ent.resize(it-ent.begin());        
//         // entity nj is connected to interface entities
//         if (ent.size()>0) {
//             on[i] = 1;     
//             k += 1;
//         }
//     }    
//     vecrecv.resize(k*3);
//     j = 0;
//     for (i=0; i<na; i++)
//         if (on[i]==1) {
//              vecrecv[0*k+j] = entrecv[0*na+i];  // store neighboring subdomain k
//              vecrecv[1*k+j] = entrecv[1*na+i];  // store local entity numbering
//              vecrecv[2*k+j] = entrecv[2*na+i];  // store global entity numbering
//              j += 1;
//         }
//     sort_columnk(vecrecv, k, 3, 0);        
//     
//     k = 0;
//     for (i=0; i<na; i++) {
//         nj = extent[i];  // global entity numbering
//         ij = ent2ind[nj]; // index of nj
//         ni = rowent2elem[ij+1] - rowent2elem[ij];
//         for (j=0; j<ni; j++) // for each elemen neighboring the entity ei        
//             elemnb[j] = colent2elem[rowent2elem[ij]+j]; // element neighboring the entity ei                
// 
//         it = set_difference(elemnb.begin(), elemnb.end(), elempart.begin(), elempart.end(), 
//                         inserter(elem, elem.begin()));    
//         elem.resize(it-elem.begin());  
//         // At least one of neighboring elements of the entity nj is not inside the overlapped subdomain i             
//         if (elem.size()>0) {
//             on[i] = 1;     
//             k += 1;
//         }                 
//         else
//             on [i] = 0;
//     }                
//     matrecv.resize(k*3);
//     j = 0;
//     for (i=0; i<na; i++)
//         if (on[i]==1) {
//              matrecv[0*k+j] = entrecv[0*na+i];  // store neighboring subdomain k
//              matrecv[1*k+j] = entrecv[1*na+i];  // store local entity numbering
//              matrecv[2*k+j] = entrecv[2*na+i];  // store global entity numbering
//              j += 1;
//         }            
//     sort_columnk(matrecv, k, 3, 0);
// }
// 
// void mkelconhdg(elcon,t2t,t2f,nproc,preconditioner,elconhdg,my_rank)
// {
//     
//     ne = elempart.size();
//     npf = param[0];
//     nfe = param[1];
//         
//     t2f.resize(ne*nfe);
//     for(j=0;j<nfe;j++) 
//         for (i=0;i<ne;i++) {
//             en = t2fhdg[j*ne+elempart[i]];
//             it = find(entpart.begin(), entpart.end(), en);
//             if (it != entpart.end())
//                 t2f[j*ne+i] = distance(entpart.begin(), it);
//             else
//                 error('Error #1 in domaindecomposition.cpp::mkelconhdg');  
//         }
//            
//     elcon.resize(npf*nfe*ne);    
//     for(i=0;i<ne;i++) 
//         for (j=0;j<nfe;j++) {
//             for (k=0;k<npf;k++) {
//                 ii[k] = elconhdg[elempart[i]*nfe*np+j*npf+k] - t2fhdg[j*ne+elempart[i]]*npf;
//                 ind[k] = t2f[j*ne+i]*npf+k;
//             }
//             for (k=0;k<npf;k++)
//                 elcon[i*nfe*npf+j*npf+k] = ind[ii[k]];
//         }
//     
//     
//     // local entity-to-element and entity-to-entity connectivities
//     elcon2entcon(bcrs_rowent2elem,bcrs_colent2elem,bcrs_rowent2ent,bcrs_colent2ent,t2f);                        
// }
// 
// void mkelconedg(elcon,t2t,t2f,nproc,preconditioner,elconhdg,my_rank)
// {
//     
//     ne = elempart.size();
//     npf = param[0];
//     nfe = param[1];
//         
//     t2f.resize(ne*nfe);
//     for(j=0;j<nfe;j++) 
//         for (i=0;i<ne;i++) 
//             t2f[j*ne+i] = t2fedg[j*ne+elempart[i]];
// 
//     vector<Int> facesInProcessor;
//     unique_copy(t2f.begin(),t2f.end(),facesInProcessor.begin());
//     
//     Int maxface = facesInProcessor.back();
//     Int nface = facesInProcessor.size();
//     vector<Int> facemapping(maxface,0);
//     for (i=0;i<nface;i++)
//         facemapping[facesInProcessor[i]] = i;
//     for(j=0;j<nfe;j++) 
//         for (i=0;i<ne;i++) 
//             t2f[j*ne+i] = facemapping[t2f[j*ne+i]];
//     
//     
//         % global element-to-entity connectivity
//         dmd{i}.elcon = elcon(:,dmd{i}.elempart); 
//         % mapping from global to local
//         entMapping = zeros(max(dmd{i}.entpart),1);
//         entMapping(dmd{i}.entpart) = 1:length(dmd{i}.entpart);
//         % local element-to-entity connectivity
//         dmd{i}.elcon = entMapping(dmd{i}.elcon);
//         if min(dmd{i}.elcon(:)) < 1; error('Something wrong.'); end
//     
//     elcon.resize(npf*nfe*ne);    
//     for(i=0;i<ne;i++) 
//         for (j=0;j<nfe*npf;j++)             
//             elcon[i*nfe*npf+j] = elconedg[elempart[i]*nfe*npf+j];      
//         
//     Int maxent = entpart.back();
//     Int nent = entpart.size();
//     vector<Int> entmapping(maxent,0);
//     for (i=0;i<nent;i++)
//         entmapping[entpart[i]] = i;
//     for(i=0;i<ne;i++) 
//         for (j=0;j<nfe*npf;j++)             
//             elcon[i*nfe*npf+j] = entmapping[elcon[i*nfe*npf+j]];
//         
//     // local entity-to-element and entity-to-entity connectivities
//     elcon2entcon(bcrs_rowent2elem,bcrs_colent2elem,bcrs_rowent2ent,bcrs_colent2ent,elcon);                        
// }
// 
// vector<vector<Int>> void getdatampi(vector<Int> & datasend, vector<Int> & nbsd)
// {
//     // send datasend to neighboring cpus 
//     // each cpu will receive and return datarecv
//     
//     Int nnbsd = nbsd.size();
//                 
//     /* blocking send */    
//     Int nsend = datasend.size();
//     for (n=0; n<nnbsd; n++) {
//         // neighboring cpu
//         neighbor = nbsd[n];
//         // send datasend to neighbor
//         MPI_Send(&datasend[0], nsend, MPI_INT, neighbor, 0, MPI_COMM_WORLD);                                
//     }
//         
//     vector<vector<Int>> datarecv;    
//     // initialize datarecv for all neighboring cpus   
//     datarecv.resize(nnbsd);
//     
//     MPI_Status status;
//     int nrecv;
//     /* blocking receive */
//     for (n=0; n<nnbsd; n++) {
//         // neighboring cpu
//         neighbor = nbsd[n];         
//         
//         // Probe for an incoming message from neighboring cpu
//         MPI_Probe(neighbor, 0, MPI_COMM_WORLD, &status);
//         
//         // When probe returns, the status object has the size and other
//         // attributes of the incoming message. Get the message size
//         MPI_Get_count(&status, MPI_INT, &nrecv);
//         
//         // Allocate a buffer to hold the incoming numbers
//         datarecv[n].resize(nrecv);
//         
//         // Now receive the message with the allocated buffer
//         MPI_Recv(&datarecv[n][0], nrecv, MPI_INT, neighbor, 0,
//                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);        
//     }      
//     
//     return datarecv;
// }
//         
// 
// void ddmpi(elist,epts,datarecv,epart,nbsd)
// {
//     
// //     vector<vector<Int>> elemrecv_fromneighbors = getdatampi(elemrecv, nbsd);
// //     vector<vector<Int>> entrecv_fromneighbors = getdatampi(entrecv, nbsd);
// //     vector<vector<Int>> vecrecv_fromneighbors = getdatampi(vecrecv, nbsd);
// //     vector<vector<Int>> matrecv_fromneighbors = getdatampi(matrecv, nbsd);
// //     vector<vector<Int>> rowent2ent_fromneighbors = getdatampi(bcrs_rowent2ent, nbsd);
// //     vector<vector<Int>> colent2ent_fromneighbors = getdatampi(bcrs_colent2ent, nbsd);        
//      
//     Int nnbsd = nbsd.size();    
//     Int nrecv, j = 0;   
//     for (n=0; n<nnbsd; n++) {
//         // neighboring cpu
//         neighbor = nbsd[n];        
//         nrecv = datarecv[n].size();
//         nrows = round(nrecv/3.0);        
//         for (i=0; i<nrows; i++) {
//             cpui = datarecv[n][0*nrows+i];
//             if (cpui==my_rank) {
//                 nlist[j] = neighbor;
//                 gi = datarecv[n][2*nrows+i];
//                 it = find(epart.begin(),epart.end(),gi);
//                 k = distance(epart.begin(),it);
//                 elist[j] = epart[k];
//                 j = j + 1;
//             }
//         }            
//     }    
//     
//     vector<Int> pts(nnbsd,0);
//     for (i=0; i<j; i++) 
//         for (n = 0; n<nnbsd; n++)
//             if (nlist[i] == nbsd[n])
//                 epts[n] += 1;            
// }
// 
// void matsendrecvedg
// {      
//     nent = matentrecv.size();
//     for (j=0; j<nent; j++) {
//         ncpu = cmatrecv[j]; // neighboring cpu
//         nj = lmatrecv[j]; // local entity nj          
//         gj = gmatrecv[j]; // global entity gj                          
//         
//         nn = bcrs_rowent2ent[nj+1] - bcrs_rowent2ent[nj];
//         for (k=0; k<nn; k++) // for each entity neighboring the entity nj        
//             entnb[k] = bcrs_colent2ent[bcrs_rowent2ent[nj]]+k]; // entities neighboring the entity nj                
//         
//         it = set_intersection(entnb.begin(), entnb.end(), lmatrecv.begin(), lmatrecv.end(), 
//                         inserter(nbentj, nbentj.begin()));    
//         nbentj.resize(it-nbentj.begin());                
//         sort(nbentj.begin(),nbentj.end());        
//         for (k=0; k<nbentj.size(); k++)
//             if (nbentj[k] == nj)
//                 nbentj.erase(nbentj.begin()+k);
//         nbentj.push_front(nj);
//         
//         for (k=0; k<nbentj.size(); k++) {
//             it = find(entnb.begin(),entnb.end(),nbentj[k]);
//             in = distance(it,entnb.begin());
//             rowj[k] = bcrs_rowent2ent[nj]+in;
//         }
//         
//         // add [ncpu rowj(1)-1] to matrecv
//         matrecv0[irecv] = ncpu;
//         matrecv1[irecv] = rowj[0];
//         irecv += 1;                        
//         
//        
//         // find neighboring elements of the entity
//         ij = ent2ind[gj];
//         nn = bcrs_rowent2elem[ij+1] - bcrs_rowent2elem[ij];
//         // list of global neighboring elements to gj
//         for (k=0; k<nn; k++) // for each elemennt neighboring the entity gj        
//             nbelemj[k] = bcrs_colent2elem[bcrs_rowent2elem[ij]]+k]; // elements neighboring the entity gj                
//         
//         // find the index of ncpu
//         for (k=0;k<nnbsd;k++)
//             if (nbsd[k] == ncpu)
//                 icpu = k;
//         
//         // find an element in matsend{ncpu}(:,3) to match gj
//         it = find(matentsend[icpu].begin()+2*nent,matentsend[icpu].end(),gj);
//         m = distance(it,matsend[icpu].begin());
//         // local entity mj on neighboring ncpu              
//         mj = matentsend[icpu][nent+m]; 
//         nn = bcrs_nbrowent2ent[icpu][mj+1] - bcrs_nbrowent2ent[icpu][mj];
//         // list of neighboring entities of the entity gj in the ncpu                
//         for (k=0; k<nn; k++) {
//             pj[k] = bcrs_nbrowent2ent[nj]]+k;
//             mbent[k] = bcrs_nbcolent2ent[icpu][pj[k]]; 
//             gbent[k] = nbentpart[icpu][mbent[k]];
//         }
//         
//         // add [i pj[0]] to matsend on neighboring ncpu                      
//         matsend0[icpu][isend[icpu]] = i;
//         matsend1[icpu][isend[icpu]] = pj[0];
//         isend[icpu] += 1;
//         
//         for (k=1; k<nbentj.size(); k++) {
//             // global entity gk on self-cpu
//             gk = entpart[nbentj[k]];        
//             // local index of gk
//             ik = ent2ind[gk]; 
//             nn = bcrs_rowent2elem[ik+1] - bcrs_rowent2elem[ik];
//             // list of global neighboring elements to gk
//             for (m=0; m<nn; m++) // for each elemennt neighboring the entity gj        
//                 nbelemk[m] = bcrs_colent2elem[bcrs_rowent2elem[ik]]+k]; // elements neighboring the entity gj                
//             // intersect nbelemj and nbelemk to get a common set
//             it = set_intersection(nbelemj.begin(), nbelemj.end(), nbelemk.begin(), nbelemk.end(), 
//                             inserter(nbelem, nbelem.begin()));    
//             nbelem.resize(it-nbelem.begin());                
//             if (IsSubset(nbelem, elempart)) {
//                 // add [ncpu rowj(k)] to matrecv
//                 matrecv0[irecv] = ncpu;
//                 matrecv1[irecv] = rowj[k];
//                 irecv += 1;                        
//                 
//                 // find an element in gbent to match gk                
//                 it = find(gbent.begin(),gbent.end(),gk);
//                 nk = distance(it,gbent.begin());
//                 // add [i pj(nk)] to matsend
//                 matsend0[icpu][isend[icpu]] = i;
//                 matsend1[icpu][isend[icpu]] = pj[nk];
//                 isend[icpu] += 1;                                
//             }            
//         }        
//     }   
// }

#endif

