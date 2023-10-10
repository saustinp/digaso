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
#include "../gmres/gmres.h"
#include "errormsg.cpp"
#include "ioutilities.cpp"
//#include "../utilities/utilities.cpp"
//#include "../utilities/sendrecvMPI.cpp"
#include "readbinaryfiles.cpp"
//#include "initializestructs.cpp"
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
    //elemstruct elems;
    //tempstruct temps;
    
    app.my_rank = world_rank;    
    mesh.my_rank = world_rank; 
    readInput(app, master, mesh, sol, string(argv[1]));    
    if (app.nproc != world_size) {
        error("MPI world size must be equal to app.nproc");
    }
    
    writeOutput(app, master, mesh, sol, string(argv[2]));
    
    dmdstruct dmd;
    createdmdstruct(dmd, mesh);    

//     string fileoutdmd = string(argv[2]) + "dmd" + NumberToString(app.my_rank) + ".bin";            
//     writedmdstruct(fileoutdmd, dmd);    
//         
//     cout<<app.my_rank<<"   "<<mesh.ne<<"   "<<mesh.nf<<"   "<<mesh.ndofuh<<endl;      
//     setupsys(sys, mesh, app);
//     cout<<app.my_rank<<"   "<<mesh.ne<<"   "<<mesh.nf<<"   "<<mesh.ndofuh<<endl;      
            
    //initializeStructs(mesh, master, app, sol, elems, sys, temps);
    
    
//     string fileapp = string(argv[1]) + "app.bin";    
//     appstruct app;    
//     readappstruct(fileapp,app);
//     app.my_rank = world_rank;   
//     
//     if (app.nproc != world_size) {
//         printiarray(app.ndims);
//         cout<<app.nproc<<"   "<<world_size<<endl;
//         error("MPI world size must be equal to app.nproc");
//     }
//     if ((app.nproc%app.nfile) != 0) {
//         error("Number of files must be divisible by number of processors");            
//     }            
//     
//     Int nmax =  *max_element(app.elemtype.begin(), app.elemtype.end());    
//     string filemaster = string(argv[1]) + "master.bin";   
//     vector< masterstruct > master(nmax+1, masterstruct());
//     readmasterstruct(filemaster, master, app);
//          
//     // read mesh structure
//     Int nprocperfile = app.nproc/app.nfile;    
//     Int filenumber = app.my_rank/nprocperfile+1; //file number
//     
//     //cout<<app.my_rank<<"   "<<app.nproc<<"   "<<app.nfile<<"   "<<nprocperfile<<"   "<<filenumber<<endl;    
//     string filemesh = string(argv[1]) + "mesh" + NumberToString(filenumber) + ".bin";                
//     meshstruct mesh;
//     readmeshstruct(filemesh,mesh,app);
//     
//     mesh.my_rank = world_rank;           
//     Int nelem = app.elemtype.size();
//     Int maxe = *max_element(app.elemtype.begin(), app.elemtype.end()); 
//     mesh.nfes.resize(maxe+1);
//     mesh.npes.resize(maxe+1);
//     mesh.npfs.resize(maxe+1);
//     mesh.ncf.resize(maxe+1);
//     mesh.ndf.resize(maxe+1);
//     mesh.perm.resize(maxe+1);    
//     mesh.cg2dg.resize(maxe+1);
//     mesh.dg2cg.resize(maxe+1);
//     Int i, j, k, e, ncf, ndf, nfe, npf, elemNode, dgNode, cgNode, alreadyExistentCGnode;           
//     for (i=0; i<nelem; i++) {        
//         e = app.elemtype[i];
//         nfe = master[e].nfe;
//         ndf = master[e].ndf;
//         mesh.ndf[e] = ndf;
//         
//         master[e].npfix.resize(nfe+1);
//         master[e].ngfix.resize(nfe+1);
//         master[e].ngfRix.resize(nfe+1);
//         master[e].npfix[0] = 0;
//         master[e].ngfix[0] = 0;
//         master[e].ngfRix[0] = 0;
//         for (j=0; j<nfe; j++) {           
//             master[e].npfix[j+1] += master[e].npfix[j] + master[e].npf[j];                  
//             master[e].ngfix[j+1] += master[e].ngfix[j] + master[e].ngf[j];                  
//             master[e].ngfRix[j+1] += master[e].ngfRix[j] + master[e].ngfR[j];                  
//         }
//         
//         mesh.nfes[e] = master[e].nfe;
//         mesh.npes[e] = master[e].npe;
//         mesh.npfs[e] = master[e].npf;     
//         
//         // Compute ncf, mesh.dg2cg and mesh.cg2dg. TODO: Introduce this as inputs in binary file
//         ncf = 0;                
//         mesh.perm[e].resize(ndf);
//         //npf = master[e].npf[0];        
//         for (j=0; j<nfe; j++) {
//             npf = master[e].npf[j];        
//             for (k = 0; k < npf; k++)
//                 mesh.perm[e][master[e].npfix[j]+k] = master[e].perm[i][k];
//         }
//         
//         mesh.cg2dg[e].resize(ndf);     // This is more than required, but that is OK
//         mesh.dg2cg[e].resize(ndf);        
//         for (i = 0; i < ndf; i++) {
//             elemNode = mesh.perm[e][i];
// 
//             alreadyExistentCGnode = 0;
//             for (j = 0; j < ncf; j++) {
//                 if (mesh.cg2dg[e][j] == elemNode) {
//                     alreadyExistentCGnode = 1;
//                     break;
//                 }
//             }
// 
//             if (alreadyExistentCGnode == 0) {
//                 cgNode = ncf;
//                 mesh.cg2dg[e][cgNode] = elemNode;      // Temporarily, cg2dg stores cg2elemNode
//                 ncf++;
//             }
//             else if (alreadyExistentCGnode == 1) {
//                 cgNode = j;
//             }
//             mesh.dg2cg[e][i] = cgNode;
//         }
//         
//         // Convert "cg2elemNode" to "cg2dg":
//         for (i = 0; i < ncf; i++) {
//             elemNode = mesh.cg2dg[e][i];
//             for (j = 0; j < ndf; j++) {
//                 if (mesh.perm[e][j] == elemNode) {
//                     dgNode = j;
//                     break;
//                 }
//             }
//             if (j == ndf) {
//                 printf("Error when computing mesh.cg2dg and mesh.dg2cg arrays\n");
//                 exit(1);
//             }
//             mesh.cg2dg[e][i] = dgNode;
//         }
//         mesh.ncf[e] = ncf;            // TODO: Include in binary files!!        
//     }        
    
//         // Compute ncf, mesh.dg2cg and mesh.cg2dg. TODO: Introduce this as inputs in binary file
//         Int ncf = 0;
//         Int nfe = master[e].nfe;
//         Int ndf = master[e].ndf;
//         Int npf = master[e].npf[0];        
//         for (i=0; i<nfe; i++)
//             for (k = 0; k < npf; k++)
//                 mesh.perm[i*npf+k] = master[e].perm[i][k];
//                     
//         mesh.cg2dg.resize(ndf);     // This is more than required, but that is OK
//         mesh.dg2cg.resize(ndf);        
//         for (i = 0; i < ndf; i++) {
//             elemNode = mesh.perm[i];
// 
//             alreadyExistentCGnode = 0;
//             for (j = 0; j < ncf; j++) {
//                 if (mesh.cg2dg[j] == elemNode) {
//                     alreadyExistentCGnode = 1;
//                     break;
//                 }
//             }
// 
//             if (alreadyExistentCGnode == 0) {
//                 cgNode = ncf;
//                 mesh.cg2dg[cgNode] = elemNode;      // Temporarily, cg2dg stores cg2elemNode
//                 ncf++;
//             }
//             else if (alreadyExistentCGnode == 1) {
//                 cgNode = j;
//             }
//             mesh.dg2cg[i] = cgNode;
//         }
//         
//         // Convert "cg2elemNode" to "cg2dg":
//         for (i = 0; i < ncf; i++) {
//             elemNode = mesh.cg2dg[i];
//             for (j = 0; j < ndf; j++) {
//                 if (mesh.perm[j] == elemNode) {
//                     dgNode = j;
//                     break;
//                 }
//             }
//             if (j == ndf) {
//                 printf("Error when computing mesh.cg2dg and mesh.dg2cg arrays\n");
//                 exit(1);
//             }
//             mesh.cg2dg[i] = dgNode;
//         }
//         mesh.ncf = ncf;            // TODO: Include in binary files!!        
            
//     // Compute mesh.isEDGelement:
//     mesh.isEDGelement.resize(mesh.ne);
//     for (i = 0; i < mesh.ne; i++) {
//         e = mesh.elementtype[i];        
//         ndf = master[e].ndf;        
//         mesh.isEDGelement[i] = 1;
//         for (j = 0; j < ndf; j++) {
//             if (mesh.elcon[ndf*i+mesh.cg2dg[e][mesh.dg2cg[e][j]]] != mesh.elcon[ndf*i+j]) {
//                 mesh.isEDGelement[i] = 0;
//                 break;
//             }
//         }
//     }        
//     
//     string filesol = string(argv[1]) + "sol" + NumberToString(filenumber) + ".bin";                
//     solstruct sol;
//     readsolstruct(filesol,sol,app);        
    
//     string fileoutapp = string(argv[2]) + "app.bin";        
//     writeappstruct(fileoutapp,app);
//     
//     string fileoutmaster = string(argv[2]) + "master.bin";        
//     writemasterstruct(fileoutmaster,master,app);
//     
//     string fileoutmesh = string(argv[2]) + "mesh" + NumberToString(filenumber) + ".bin";                
//     writemeshstruct(fileoutmesh,mesh,app);
//     
//     string fileoutsol = string(argv[2]) + "sol" + NumberToString(filenumber) + ".bin";            
//     writesolstruct(fileoutsol,sol,app);
                               
//     dmdstruct dmd;
//     createdmdstruct(dmd, mesh);    
//     
//     string fileoutdmd = string(argv[2]) + "dmd" + NumberToString(app.my_rank) + ".bin";            
//     writedmdstruct(fileoutdmd, dmd);
    
    // Finalize the MPI environment:
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
        
    return 0;             
}
