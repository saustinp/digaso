#ifndef __MKMESH_H
#define __MKMESH_H

#include "../utilities/wrapper.h"
//#include "../utilities/datastructs.h"

// element info struct
struct elementinfo {    
    Int dim;      // spatial dimension
    Int elemtype; // type     
    
    Int nfe; // number of faces    
    Int nle;  // number of edges 
    Int nve; // number of vertices 
    Int nqf; // number of quad faces
    Int ntf; // number of tet faces
    vector< Int > facetype; // 0 triangular, 1 quadrilateral 
    vector< vector< Int > > face; // face connectivity     
    vector< vector< Int > > edge; // edge connectivity     
};

// element struct
struct elementstruct {        
    //elementinfo eleminfo;
    Int dim;      // spatial dimension
    Int elemtype; // type         
    Int nfe; // number of faces    
    Int nle;  // number of edges 
    Int nve; // number of vertices 
    vector< Int > facetype; // 0 triangular, 1 quadrilateral 
    vector< vector< Int > > face; // face connectivity     
    //vector< vector< Int > > edge; // edge connectivity     
    
    vector < Int > porder; // polynomial degrees along each direction
    Int nodetype; // 0 uniform, 1 optimal        
    Int npe; // number of nodes
    Int nte; // number of subelements
    vector< double>  plocvl; // nodal points
    vector< Int > tlocvl;    // subelements
    vector< Int > npf; // number of nodes on faces
    vector< Int > ntf; // number of subelements on faces
    vector< vector< double> > plocfc; // nodal points on faces
    vector< vector< Int > > tlocfc;   // subelements on faces
    vector< vector< Int > > perm;     // indices of the nodes on faces
};

struct masterstruct {
    vector<Int> pgauss;    
    Int ngv;
    vector<Int> ngf;    
        
    vector<double> gpvl;
    vector<double> gwvl;    
    vector<double> shapvl;
    vector<double> shapvt;
    vector<double> shapvg;
    vector<double> shapvgdotshapvl;
    vector< vector<double> > gpfc;
    vector< vector<double> > gwfc;    
    vector< vector<double> > shapfc;
    vector< vector<double> > shapft;
    vector< vector<double> > shapfg;
    vector< vector<double> > shapfgdotshapfc;
    
    vector<double> shapnv;
    vector<double> shapnvt;   
    vector<double> projLowP;
    
    // NEW FIELDS FOR NEW MATRIX ASSEMBLY:
    Int pgaussR;
    Int pgaussJ;
    Int pgaussQ;
    
    Int quadTypeR;
    Int quadTypeJ;
    Int quadTypeQ;
    
    Int nqvR;
    Int nqvQ;
    Int nqvJ;
    vector<double> gpvlR;
    vector<double> gpvlJ;
    vector<double> gpvlQ;
    vector<double> gwvlR;
    vector<double> gwvlJ;
    vector<double> gwvlQ;
    vector<double> shapvlR;
    vector<double> shapvlJ;
    vector<double> shapvlQ;
    vector<double> shapvtR;
    vector<double> shapvtJ;
    vector<double> shapvtQ;
    vector<double> shapvgR;
    vector<double> shapvgJ;
    vector<double> shapvgQ;
    vector<double> shapvgdotshapvlR;
    vector<double> shapvgdotshapvlJ;
    vector<double> shapvgdotshapvlQ;
    
    vector<Int> nqfR;
    vector<Int> nqfQ;
    vector<Int> nqfJ;
    vector< vector<double> > gpfcR;
    vector< vector<double> > gpfcJ;
    vector< vector<double> > gpfcQ;
    vector< vector<double> > gwfcR;
    vector< vector<double> > gwfcJ;
    vector< vector<double> > gwfcQ;
    vector< vector<double> > shapfcR;
    vector< vector<double> > shapfcJ;
    vector< vector<double> > shapfcQ;
    vector< vector<double> > shapftR;
    vector< vector<double> > shapftJ;
    vector< vector<double> > shapftQ;
    vector< vector<double> > shapfgR;
    vector< vector<double> > shapfgJ;
    vector< vector<double> > shapfgQ;
    vector< vector<double> > shapfgdotshapfcR;
    vector< vector<double> > shapfgdotshapfcJ;
    vector< vector<double> > shapfgdotshapfcQ;
};
       
struct meshstruct {
    Int dim;   // spatial dimension    
    Int ne; // number of elements
    Int nf; // number of faces
    Int nv; // number of vertices   
    Int nfemax; // maximum number of faces on elements    
    Int nlemax;  // maximum number of edges on elements
    Int nvemax; // maximum number of vertices on elements
    Int nvfmax; // maximum number of vertices on faces
    Int npemax; // maximum number of nodes on elements
    Int npfmax; // maximum number of nodes on faces    
    Int nodetype; // 0 uniform, 1 optimal        
    Int hybrid;
            
    vector<double> p; // coordinates of vertice
    vector<double> dgnodes; // coordinates of DG nodes   
    vector<Int> t;  // element-to-vertice connectivities
    vector<Int> dgnodesidx; // element index of dgnodes        
    vector<Int> elcon; // element-to-entity connectivities
    vector<Int> elconidx; // element index of elcon        
    vector< vector<Int> > npfidx;  // face index of npf
    vector<Int> bf;  // numbering of boundary faces
    vector<Int> t2f; // element-to-face connectivities
    vector<Int> t2t; // element-to-element connectivities
    vector<Int> f;   // face-to-element connnectivities 
    vector<Int> cg2dg;
    vector<Int> dg2cg;    
     
    vector<Int> porder;    
    vector<Int> pgauss;
    vector<Int> quadType;    
    vector<elementstruct> elements; // types of elements in the mesh    
    vector<masterstruct> masters; // types of masters in the mesh    
    vector<Int> elementtype;  // associated with the structure elements 
    vector<Int> physics;     // associated with the governing equations            
    vector<Int> isEDGface;    
    
    vector<double> elemMeasure; // Element measure: ne
    vector<double> hAvg;        // Characteristic element size: ne
    vector<double> M;           // Metric tensor at dgnodes
    vector<double> Minv;        // Inverse of metric tensor at dgnodes
};

void mkmesh(meshstruct &mesh, vector<Int> &porder, Int dim, Int nodetype);

#endif
