#ifndef __DMDLIB
#define __DMDLIB

// Written by: C. Nguyen & P. Fernandez

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <mpi.h>

#ifndef HAVE_MPI
#define HAVE_MPI
#endif

using namespace std;

#include "digasopre.h"
#include "errormsg.cpp"
#include "mkdmd.cpp"

void createdmdstruct(dmdstruct &dmd, meshstruct &mesh)
{
    Int N = mesh.ndims.size();
    Int nie = mesh.ndims[17];  // number of interior elements
    Int nif = mesh.ndims[18];  // number of interior faces                    
    Int nin = mesh.ndims[19];  // number of interior entities                     
    Int nxe = mesh.ndims[21];  // number of exterior elements
    Int nxf = mesh.ndims[22];  // number of exterior faces     
    Int nxn = mesh.ndims[23];  // number of exterior entities            
    Int hybrid = mesh.ndims[N-2];        
    
    if (hybrid == 0) {  // HDG              
        domaindecomposition(dmd, mesh.t2f, mesh.elcon, mesh.t2f, mesh.elementtype, mesh.extintelem, mesh.extintface,
                mesh.elem2cpu, mesh.face2cpu, mesh.nfes, mesh.npfs, nie, nin, nxe, nxn, 1, mesh.my_rank);                        
        Int nnbsd = (Int) dmd.nbsd.size();                  
        vector< dmdstruct > nbdmd(nnbsd, dmdstruct());
        make_nbdmd(nbdmd, dmd);        
        sendrecvhdg(dmd);
    }
    else {                        
        domaindecomposition(dmd, mesh.elcon, mesh.elcon, mesh.t2f, mesh.elementtype, mesh.extintelem, mesh.extintent,
                mesh.elem2cpu, mesh.ent2cpu, mesh.nfes, mesh.npfs, nie, nin, nxe, nxn, 0, mesh.my_rank);                
        Int nnbsd = (Int) dmd.nbsd.size();          
        vector< dmdstruct > nbdmd(nnbsd, dmdstruct());
        make_nbdmd(nbdmd, dmd);        
        sendrecvedg(dmd, nbdmd);
    }                            
}    

#endif
