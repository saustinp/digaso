#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <stdlib.h>

// Written by: C. Nguyen & P. Fernandez

using namespace std;

//#include "domaindecomposition.h"
#include "errormsg.cpp"
#include "domaindecomposition.cpp"
#include "ioutil.cpp"

int main(int argc, char** argv) 
{
   
    if( argc == 3 ) {
    }
    else {
      printf("Usage: ./cppfile InputFile OutputFile\n");
      return 1;
    }            
    
    string filein = string(argv[1]) + ".bin";       // Name of binary file with input data
    string fileout = string(argv[2]) + ".bin";        // Name of binary file to write the solution
    
    vector<Int> param, elcon, t2t, t2f, elemall, entall;        
    readInput(filein, param, elcon, t2t, t2f, elemall, entall);
    
    /* Get dimensions from param array */
    Int ne = param[0];
    Int nf = param[1];
    Int ndof = param[2];
    Int nfe = param[3];
    Int nee = param[4];
    Int npe = param[5];
    Int npf = param[6];                
    Int preconditioner = param[10];
    Int overlappinglevel = param[11]; 
    Int hdg = param[12];
    Int nproc = param[13];
    Int np = npf*nfe;
    
    Int i, j, k, m, n, my_rank;
    vector<Int> elem2ent(ne*nfe);        
    
    if (hdg == 1) {
        // take transpose        
        for (i=0;i<ne;i++)
            for(j=0;j<nfe;j++)         
                elem2ent[i*nfe+j] = t2f[j*ne+i];        
    }
    else
        elem2ent = elcon;
    
//     print1iarray(&param[0], 20);
//     //print2iarray(&elcon[0], npf*nfe, ne);    
//     print2iarray(&t2t[0], ne, nee);
//     print2iarray(&t2f[0], ne, nfe);
//     print1iarray(&elemall[0], (Int) elemall.size());
//     print1iarray(&entall[0], (Int) entall.size());
    
    vector<Int> extintelem, intelem, intent, intfent, elemintfent, intfelem, extelem;
    vector<Int> rowent2elem, colent2elem, rowent2ent, colent2ent, ent2ind;
    vector<Int> elempart, elempartpts, extent, entpart, entpartpts, elemrecv;
    vector<Int> entrecv, nbsd, vecrecv, matrecv, dmdelcon, dmdt2f;    
    vector<Int> bcrs_rowent2elem, bcrs_colent2elem, bcrs_rowent2ent, bcrs_colent2ent;
    vector<dmdstruct> dmdall(nproc, dmdstruct());
    
    // list of  elements in nonoverlapping subdomain my_rank   
    for (i = 0; i<nproc; i++) {
        my_rank = i;
        
//         intelem = find(elemall,my_rank); 
//     
//         // list of  entities in nonoverlapping subdomain my_rank                       
//         intent = find(entall,my_rank);
    
        cout<< "rank: " << my_rank << endl;
//         print1iarray(&intelem[0], (Int) intelem.size());
//         print1iarray(&intent[0], (Int) intent.size());
        
        // list of elements on overlapping subdomain my_rank        
        //overlappingsubdomain(extintelem,intelem,t2t,overlappinglevel,ne,nee);   
        
        //vector<Int> rowent2elem, colent2elem, rowent2ent, colent2ent, ent2ind;
        entelempart(intelem, intent, rowent2elem, colent2elem, rowent2ent, colent2ent, 
                    ent2ind, elem2ent, t2t, elemall, entall, np, ne, nee, overlappinglevel, my_rank);                
        
//         cout<<rowent2elem.size()<<endl;
//         cout<<colent2elem.size()<<endl;
//         cout<<rowent2ent.size()<<endl;
//         cout<<colent2ent.size()<<endl;
//         cout<<ent2ind.size()<<endl;
        
        find_intfent(intfent, elemintfent, intent, ent2ind, rowent2elem, colent2elem, rowent2ent, colent2ent);
        
//         print1iarray(&intfent[0], (Int) intfent.size());
//         print1iarray(&elemintfent[0], (Int) elemintfent.size());

        find_extelem(intfelem, extelem, elemintfent, intelem, t2t, preconditioner, ne, nee);
     
        //print1iarray(&intfelem[0], (Int) intfelem.size());
        //print1iarray(&extelem[0], (Int) extelem.size());

        make_elempart(elempart, elempartpts, intelem, intfelem, extelem);      
        
        //printiarray(elempart);
        //printiarray(elempartpts);
        
        find_extent(extent, intfelem, extelem, intent, elem2ent, np);
        
        //printiarray(extent);        
        //cout<<extent.size()<<endl;
        
        make_elempart(entpart, entpartpts, intent, intfent, extent);
        
        //printiarray(entpart);
        //printiarray(entpartpts);
                 
        find_elemrecv(elemrecv, extelem, intelem, elemall);
        
        //n = round(((double) elemrecv.size())/3);
        //print2iarray(&elemrecv[0], n, 3);
        
        find_entrecv(entrecv, nbsd, extent, intent, entall);
        
        //n = round(((double) entrecv.size())/3);
        //print2iarray(&entrecv[0], n, 3);
        //printiarray(nbsd);
        
        find_vecrecv(vecrecv, extent, intfent, entrecv, ent2ind, rowent2ent, colent2ent);
        
        //n = round(((double) vecrecv.size())/3);
        //print2iarray(&vecrecv[0], n, 3);
        
        find_matrecv(matrecv, extent, elempart, entrecv, ent2ind, rowent2elem, colent2elem);   
            
//         n = round(((double) matrecv.size())/3);
//         print2iarray(&matrecv[0], n, 3);
//         cout<<n<<"  ,  "<<3<<endl;
        
        if (hdg==1)
            make_elconhdg(dmdelcon, dmdt2f, bcrs_rowent2elem, bcrs_colent2elem, 
                bcrs_rowent2ent, bcrs_colent2ent, entpart, elempart, 
                elcon, t2f, nfe, npf);
        else
            make_elconedg(dmdelcon, dmdt2f, bcrs_rowent2elem, bcrs_colent2elem, 
                bcrs_rowent2ent, bcrs_colent2ent, entpart, elempart, 
                elem2ent, t2f, nfe, npf);

        //printiarray(bcrs_rowent2elem);
        //printiarray(bcrs_colent2elem);
        printiarray(bcrs_rowent2ent);
        //printiarray(bcrs_colent2ent);
        
        domaindecomposition(dmdall[i],elem2ent,elcon, t2t, t2f, elemall, entall, 
                preconditioner, overlappinglevel, hdg, my_rank);
        
        printiarray(dmdall[i].elempart);
    }
    
//     vector<Int> extintelem, intelem;
//     overlappingsubdomain(extintelem,intelem,t2t,overlappinglevel,ne,nee);
//     
//     ofstream out(fileout.c_str(), ios::out | ios::binary);
//     if (!out)
//         exit(-1);
//     writeiarray(out, extintelem, (Int) extintelem.size());
//     writeiarray(out, intelem, (Int) intelem.size());
//     //writedarray(out, master.shapfcR, (Int) master.shapfcR.size());
//     out.close();    
// 
    return 0;     
    
// void domaindecomposition(dmdstruct &dmd, vector<Int> &elcon, vector<Int> &elconhdg, 
//         vector<Int> &t2t, vector<Int> &t2f, vector<Int> &elemall, vector<Int> &entall, 
//         Int preconditioner, Int overlappinglevel, Int hdg, Int my_rank)
    
}
