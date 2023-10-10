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
//#include "../gmres/gmres.h"
#include "preprocessing.cpp"
#include "../assembly/assembly.cpp"

//#include "../utilities/utilities.cpp"
//#include "../utilities/sendrecvMPI.cpp"
//#include "mkdmd.cpp"
//#include "mkmesh.cpp"
// #include "../assembly/elementassemblystandard.cpp"
// #include "../assembly/elementschur.cpp"
// #include "../assembly/globalassembly.cpp"
//#include "../assembly/assembly.cpp"
//#include "../solvers/dgsolversmpi.cpp"
//#include "../gmres/gmresmpilib.cpp"

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
    
    appstruct app;
    vector< masterstruct > master;
    meshstruct mesh;
    solstruct sol;
    sysstruct sys;
    elemstruct elems;
    tempstruct temps;
    
    app.debugmode = 1;
    app.my_rank = world_rank;    
    mesh.my_rank = world_rank; 
    sys.my_rank = world_rank; 
    sys.noThreads = 1;
    preprocessing(mesh, master, app, sol, elems, sys, temps, string(argv[1]));    
    
    double elemtimes[32];
    Int ie = 0;
    Int e = mesh.elementtype[ie];
    assembleElementMatrixVector(elems, mesh, master[e], app, sol, temps, ie, &elemtimes[0]);
    
    // Finalize the MPI environment:
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
        
    return 0;             
}
