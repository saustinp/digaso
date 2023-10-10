#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

// Written by: C. Nguyen & P. Fernandez

#ifndef HAVE_MPI
#define HAVE_MPI
#endif

using namespace std;

#include "digasopre.h"
#include "dmd.h"
#include "errormsg.cpp"
#include "ioutilities.cpp"
#include "readbinaryfiles.cpp"
//#include "mkdmd.cpp"
//#include "mkmesh.cpp"

int main(int argc, char** argv) 
{
   
    if( argc == 3 ) {
    }
    else {
      printf("Usage: ./cppfile InputFile OutputFile\n");
      return 1;
    }            
    
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    string fileapp = string(argv[1]) + "app.bin";    
    appstruct app;    
    readappstruct(fileapp,app);
    app.my_rank = world_rank;   
    
    if (app.nproc != world_size) {
        printiarray(app.ndims);
        cout<<app.nproc<<"   "<<world_size<<endl;
        error("MPI world size must be equal to app.nproc");
    }
    if ((app.nproc%app.nfile) != 0) {
        error("Number of files must be divisible by number of processors");            
    }            
    
    Int nmax =  *max_element(app.elemtype.begin(), app.elemtype.end());    
    string filemaster = string(argv[1]) + "master.bin";   
    vector< masterstruct > master(nmax+1, masterstruct());
    readmasterstruct(filemaster, master, app);
     
    // read mesh structure
    Int nprocperfile = app.nproc/app.nfile;    
    Int filenumber = app.my_rank/nprocperfile+1; //file number
    
    //cout<<app.my_rank<<"   "<<app.nproc<<"   "<<app.nfile<<"   "<<nprocperfile<<"   "<<filenumber<<endl;    
    string filemesh = string(argv[1]) + "mesh" + NumberToString(filenumber) + ".bin";                
    meshstruct mesh;
    readmeshstruct(filemesh,mesh,app);        
    mesh.my_rank = world_rank;       
    
    Int nelem = app.elemtype.size();
    Int maxe = *max_element(app.elemtype.begin(), app.elemtype.end()); 
    mesh.nfes.resize(maxe+1);
    mesh.npes.resize(maxe+1);
    mesh.npfs.resize(maxe+1);
    Int i, e;
    for (i=0; i<nelem; i++) {
        e = app.elemtype[i];
        mesh.nfes[e] = master[e].nfe;
        mesh.npes[e] = master[e].npe;
        mesh.npfs[e] = master[e].npf;
    }        
        
    string filesol = string(argv[1]) + "sol" + NumberToString(filenumber) + ".bin";                
    solstruct sol;
    readsolstruct(filesol,sol,app);        
                
    dmdstruct dmd;
    createdmdstruct(dmd, mesh);    
    
    string fileoutapp = string(argv[2]) + "app.bin";        
    writeappstruct(fileoutapp,app);
    
    string fileoutmaster = string(argv[2]) + "master.bin";        
    writemasterstruct(fileoutmaster,master,app);
    
    string fileoutmesh = string(argv[2]) + "mesh" + NumberToString(filenumber) + ".bin";                
    writemeshstruct(fileoutmesh,mesh,app);
    
    string fileoutsol = string(argv[2]) + "sol" + NumberToString(filenumber) + ".bin";            
    writesolstruct(fileoutsol,sol,app);
                
    string fileoutdmd = string(argv[2]) + "dmd" + NumberToString(app.my_rank) + ".bin";            
    writedmdstruct(fileoutdmd, dmd);
    
    // Finalize the MPI environment:
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
        
    return 0;             
}
