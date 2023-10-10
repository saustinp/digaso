#ifndef __MKMESH
#define __MKMESH
        
#include "mkmesh.h"
#include "mathutilities.cpp"
#include "mkmasternodes.cpp"
#include "mkelement.cpp"
#include "mkdgnodes.cpp"
#include "mkt2f.cpp"
#include "mkmaster.cpp"
#include "mkelcon.cpp"

// p t elementtype dim // porder pgauss nodetype quadType hybrid
void mkmesh(meshstruct &mesh, vector<Int> &porder, Int nodetype, vector<Int> &pgauss, vector<Int> &quadType, Int hybrid)
{       
    Int dim = mesh.dim;
    mesh.porder = porder;
    mesh.nodetype = nodetype;
    mesh.pgauss = pgauss;
    mesh.quadType = quadType;
    mesh.hybrid = hybrid;
    
    Int i, j, n, ei, ne, nfe, nmax = 5;    
    vector<elementstruct> elements(nmax, elementstruct());    
    vector< vector<double> > philocvl(nmax, vector<double>());   
    vector< vector< vector<double> > > philocfc(nmax, vector< vector<double> >());       
    vector <vector<Int> > npfidx(nmax, vector<Int>());     
    
    //
    mesh.masters.resize(nmax);  
    vector<Int> elem = uniqueiarray(mesh.elementtype);            
    Int nelem = elem.size();                            
    for (n=0; n<nelem; n++) {
        i = elem[n];
        
        elements[i] = mkelementstruct(dim, i, porder, nodetype);        
        localbasiselemface(philocvl[i], philocfc[i], elements[i].plocvl, elements[i].plocfc, dim, elements[i].elemtype, elements[i].nve);                                    
        mkmaster(mesh.masters[i], elements[i], pgauss, quadType);        
        
        nfe = elements[i].nfe;
        npfidx[i].resize(nfe+1);
        npfidx[i][0] = 0;
        for (j=0; j<nfe; j++)
            npfidx[i][j+1] = npfidx[i][j] + elements[i].npf[j];
    }     
    mesh.elements = elements;
    mesh.npfidx = npfidx;
        
    //
    mesh.nvfmax = dim;
    Int maxelem = *max_element(elem.begin(),elem.end());
    if ((maxelem>0) & (dim==3))
       mesh.nvfmax = 4;          
    mesh.nvemax = 0;
    mesh.nfemax = 0;    
    mesh.npemax = 0;
    mesh.npfmax = 0;    
    for (n=0; n<nelem; n++) {
        i = elem[n];
        if (elements[i].npe > mesh.npemax)
            mesh.npemax = elements[i].npe;
        if (elements[i].nve > mesh.nvemax)
            mesh.nvemax = elements[i].nve;
        if (elements[i].nfe > mesh.nfemax)
            mesh.nfemax = elements[i].nfe;
        maxelem = *max_element(elements[i].npf.begin(),elements[i].npf.end());
        if (maxelem > mesh.npfmax)
            mesh.npfmax = maxelem;
    }    
        
    // TO DO: parallel version of this function
    mkt2f(mesh.t2f, mesh.t2t, mesh.f, mesh.t, mesh.elementtype, dim);  
    nmax = mesh.nvfmax+3;
    Int nf = round(((double) mesh.f.size())/(nmax));  
    mesh.isEDGface.resize(nf,0);
    if ((hybrid==0) | (hybrid==1)) { // HDG or EDG
        for (i=0; i<nf; i++)
            mesh.isEDGface[i] = hybrid;
    }
    else { // IEDG
        for (i=0; i<nf; i++) {
            if (mesh.f[i*nmax+nmax-1] < 0)
                mesh.isEDGface[i] = 0;  // HDG faces on the domain boundary
            else
                mesh.isEDGface[i] = 1;  // HDG faces on the domain boundary
        }
    }
        
    // TO DO: Parallel version of this function
    createnodes(mesh.dgnodes, mesh.dgnodesidx, mesh.p, mesh.t, mesh.elementtype, elements, philocvl, dim);                  
    
    // TO DO: Parallel version of this function
    ne = mesh.elementtype.size();        
    mesh.elconidx.resize(ne+1,0);
    mesh.elconidx[0] = 0;
    for (i=0; i<ne; i++) {
        ei = mesh.elementtype[i];        
        mesh.elconidx[i+1] = mesh.elconidx[i] + accumulate(elements[ei].npf.begin(), elements[ei].npf.end(), 0);
    }        
    mesh.elcon.resize(mesh.elconidx[ne]);        
    mkelcon(mesh.elcon, mesh.elconidx, mesh.npfidx, mesh.dgnodes, mesh.dgnodesidx, 
            mesh.elements, mesh.p, philocfc, mesh.f, mesh.t2f, mesh.elementtype, mesh.isEDGface,
            mesh.dim, mesh.nvfmax, mesh.npfmax, mesh.nfemax);      
        
    mesh.ne = ne;
    mesh.nf = mesh.isEDGface.size();
}
       
// struct meshstruct {
//     Int dim;   // spatial dimension    
//     Int ne; // number of elements
//     Int nf; // number of faces
//     Int nv; // number of vertices   
//     Int nfemax; // maximum number of faces on elements    
//     Int nlemax  // maximum number of edges on elements
//     Int nvemax; // maximum number of vertices on elements
//     Int nvfmax; // maximum number of vertices on faces
//     Int npemax; // maximum number of nodes on elements
//     Int npfmax; // maximum number of nodes on faces    
//     Int nodetype; // 0 uniform, 1 optimal        
//     
//     vector<double> p; // coordinates of vertice
//     vector<double> dgnodes; // coordinates of DG nodes   
//     vector<Int> t;  // element-to-vertice connectivities
//     vector<Int> dgnodesidx; // element index of dgnodes        
//     vector<Int> elcon; // element-to-entity connectivities
//     vector<Int> elconidx; // element index of elcon        
//     vector< vector<Int> > npfidx;  // face index of npf
//     vector<Int> bf;  // numbering of boundary faces
//     vector<Int> t2f; // element-to-face connectivities
//     vector<Int> t2t; // element-to-element connectivities
//     vector<Int> f;   // face-to-element connnectivities 
//     vector<Int> cg2dg;
//     vector<Int> dg2cg;    
// 
//     vector<Int> porder;
//     vector<elementstruct> elements; // types of elements in the mesh    
//     vector<masterstruct> masters; // types of masters in the mesh    
//     vector<Int> elementype;  // associated with the structure elements 
//     vector<Int> physics;     // associated with the governing equations            
//     vector<Int> isEDGelement;    
//     
//     vector<double> elemMeasure; // Element measure: ne
//     vector<double> hAvg;        // Characteristic element size: ne
//     vector<double> M;           // Metric tensor at dgnodes
//     vector<double> Minv;        // Inverse of metric tensor at dgnodes
// };

#endif