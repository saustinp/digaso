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


//#include "mkmesh.h"
#include "errormsg.cpp"
#include "mkmesh.cpp"
// #include "mathutilities.cpp"
// #include "mkmasternodes.cpp"
// #include "mkelement.cpp"
// #include "mkdgnodes.cpp"
// #include "mkt2f.cpp"
// #include "mkmaster.cpp"
// #include "mkelcon.cpp"
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
        
    //vector<Int> param, t, elementtype;        
    //vector<double> p;
    //readInput(filein, param, t, elementtype, p);
      
    vector<Int> param;
    meshstruct mesh;
    readInput(filein, param, mesh.t, mesh.elementtype, mesh.p);
    
    Int ne = param[0];
    Int np = param[1];        
    Int dim = param[2];        
    Int nvemax = param[3];                
    Int i;
    
    vector<Int> porder(dim,0);     
    vector<Int> pgauss(3,0);
    vector<Int> quadType(3,0);
    Int nodetype = 0;    
    Int hybrid = 0;    
    for (i=0; i<dim; i++) 
        porder[i] = 3;
    for (i=0; i<3; i++) {
        pgauss[i] = 2*porder[i];
        quadType[i] = 0;
    }
    
    mkmesh(mesh, porder, nodetype, pgauss, quadType, hybrid);
    
    return 0;         
}
