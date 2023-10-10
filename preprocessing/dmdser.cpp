#ifndef __DMDSER
#define __DMDSER

// static void uniqueiarray(vector<Int> & a, vector<Int> & b)
// {	
//     // input: integer array b
//     // output: sorted integer array a without duplicate elements
//     
//     // make a new copy of b
//     a = b;
//             
// 	// First Sort the given range to bring duplicate
// 	// elements at consecutive positions
// 	sort(a.begin(), a.end());
//   
// 	vector<Int>::iterator newEnd;
//  
// 	// Override duplicate elements
// 	newEnd = unique(a.begin(), a.end());
//  
//     // remove duplicate elements
// 	a.erase(newEnd, a.end());
// }
// 
// static void uniqueiarray(vector<Int> & a, vector<Int> & b, Int m, Int n)
// {	
//     // input: integer array b
//     // output: sorted integer array a without duplicate elements
//     
//     // make a new copy of b    
//     a.resize(n);    
//     for (Int i=m; i<m+n; i++)
//         a[i-m] = b[i];    
//             
// 	// First Sort the given range to bring duplicate
// 	// elements at consecutive positions
// 	sort(a.begin(), a.end());
//   
// 	vector<Int>::iterator newEnd;
//  
// 	// Override duplicate elements
// 	newEnd = unique(a.begin(), a.end());
//  
//     // remove duplicate elements
// 	a.erase(newEnd, a.end());
// }

static void elcon2entcon(vector<Int> &rowent2elem,vector<Int> &colent2elem,vector<Int> &rowent2ent,vector<Int> &colent2ent,
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

void dmdserial(sysstruct &sys, meshstruct &mesh, appstruct &app)
{
    Int N = mesh.ndims.size();
    Int hybrid = mesh.ndims[N-2];        
    
    mesh.ne = mesh.elementtype.size();    
    mesh.nf = *max_element(mesh.t2f.begin(), mesh.t2f.end())+1;
    mesh.ndofuh = *max_element(mesh.elcon.begin(), mesh.elcon.end())+1;
    
    vector<Int> rowent2elem, colent2elem;
    if (hybrid == 0) { // HDG                
        elcon2entcon(rowent2elem,colent2elem,sys.ent2entStart,sys.ent2ent,mesh.t2f,mesh.elementtype,mesh.nfes,mesh.nfemax,mesh.ne);                 
    }
    else {
        Int nmax = mesh.npfmax*mesh.nfemax;    
        vector<Int> nmaxs = mesh.nfes;
        Int nfe, nelem = mesh.nfes.size();
        for(Int i=0;i<nelem;i++) {
            nfe = mesh.nfes[i];
            nmaxs[i] = 0;
            for (Int j=0; j<nfe; j++)
                nmaxs[i] += mesh.npfs[i][j];                 
        }
        elcon2entcon(rowent2elem,colent2elem,sys.ent2entStart,sys.ent2ent,mesh.elcon,mesh.elementtype,nmaxs,nmax,mesh.ne);         
    }
    
    sys.nproc = app.nproc;
    if (hybrid == 0) //HDG
        sys.blkSize = app.nch*mesh.npfmax;
    else
        sys.blkSize = app.nch;    
    sys.BJ_nrows = sys.ent2entStart.size()-1;
    sys.numEntities = sys.ent2entStart.size()-1;
    sys.numBlocks = sys.ent2ent.size();    
        
    Int a;
    sys.maxBlocksPerRow = 0;
    sys.minBlocksPerRow = 9999999;
    for (Int i=1; i<sys.ent2entStart.size(); i++) {
        a = sys.ent2entStart[i]-sys.ent2entStart[i-1];
        if (a > sys.maxBlocksPerRow)
            sys.maxBlocksPerRow = a;
        if (a < sys.minBlocksPerRow)
            sys.minBlocksPerRow = a;
    }    
}

#endif
