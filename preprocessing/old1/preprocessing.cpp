#ifndef __PREPROCESSING
#define __PREPROCESSING

#include "readbinaryfiles.cpp"
#include "initializestructs.cpp"

void preprocessing(meshstruct &mesh, vector< masterstruct > &master, appstruct &app, solstruct &sol, elemstruct &elems, sysstruct &sys, tempstruct &temps, string filein) 
{   
    // read from binary files to construct app, master, mesh, sol
    readInput(app, master, mesh, sol, filein);                    
    
    // setup sys struct
    setupsys(sys, mesh, app);    
            
    // allocate memory 
    initializeStructs(mesh, master, app, sol, elems, sys, temps);        
}

#endif