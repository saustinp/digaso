#ifndef __DIGASOPRE_H
#define __DIGASOPRE_H

#include "../utilities/wrapper.h"

//appstruct
struct appstruct {
    Int my_rank;
    Int nproc;
    Int nfile;
    vector<Int> ndims;
    vector<Int> meshsize;
    vector<Int> solsize;
    vector<Int> porder;
    vector<Int> elemtype;
    vector<Int> nodetype;
    vector<Int> pgauss;
    vector<Int> pgaussR;
    vector<Int> quadtype;
    vector<Int> bcm;
    vector<double> bcs;
    vector<Int> bcd;
    vector<double> bcv;
    vector<double> dt;      
    vector<Int> flag;
    vector<double> factor;
    vector<Int> problem;
    vector<double> physicsparam;
    vector<double> solversparam;
    
    double time;            /* current time */
    double fc_u;            /* factor when discretizing the time derivative of the U equation. Allow scalar field for local time stepping in steady problems? */
    double fc_q;            /* factor when discretizing the time derivative of the Q equation. Allow scalar field for local time stepping in steady problems? */
    double fc_p;            /* factor when discretizing the time derivative of the P equation. Allow scalar field for local time stepping in steady problems? */
    double dtfc;            /* needed only for augmented Lagrangian */
    double alpha;           /* needed only for augmented Lagrangian */

    Int tdep;               /* flag to determine unsteady or steady problems */
    Int wave;               /* flag to determine wave or non-wave problems */
    Int alag;               /* flag to determine augmented lagrangian */
    Int ALEflag;            // Flag for Arbitrary Lagrangian-Eulerian (ALE) formulation. 0: No ALE; 1: Translation; 2: Translation + rotation; 3: Translation + rotation + deformation
    Int adjoint;            /* flag to determine adjoint or primal problems */
    Int linearProblem;      /* flag to determine linear or nonlinear problems */
    Int flag_q;             /* flag to determine if the Q equation is also discretized  */
    Int flag_p;             /* flag to determine if the P equation is also discretized  */
    Int flag_g;             /* flag to determine if the GCL equation is also discretized  */

    Int quasiNewton;        // flag to determine quasi Newton or full Newton
    Int reuseJacobian;      // flag to determine whether or not the Jacobian matrix needs to be computed
    Int reuseResidual;      // flag to determine whether or not the residual vector needs to be computed
    Int reusePreconditioner;// flag to determine whether or not the preconditioner needs to be computed
    Int reuseOrdering;      // flag to determine whether or not the DOF ordering needs to be computed

    Int hybrid;             // Discretization scheme. 0: HDG; 1: EDG; 2: IEDG
    Int appname;            // 0: Compressible Euler; 1: Compressible Navier-Stokes; 2: Compressible regularized RANS-SA
    Int linearSolver;       // 0: direct solver; 1: gmres; etc. (future: CG, CGS, QMR)
    Int dirkStage;          // DIRK stages
    Int dirkOrder;          // DIRK order
    Int BDFsteps;           // Number of steps for BDF scheme
    Int jacobianStep;       // reuse Jacobian for every jacobianStep time steps
    Int orderingStep;       // reuse Ordering for every orderingStep time steps

    //vector<Int> DIRKnotConverged;
    //vector<Int> BDFnotConverged;
    
    Int AVflag;             // Flag for artificial viscosity. 0: No artificial viscosity; 1: Homogeneous artificial viscosity (C. Nguyen's formulation); 2: Hypersonic homogeneous artificial viscosity (C. Nguyen's formulation)
                                    // 3: Isotropic artificial viscosity (D. Moro's formulation). 4: Latest version of the model (taking the best of all previous models)
                                    // 8: Density smoothness sensor (Per's approach)
    double rampFactor;      // Ramp factor for artificial viscosity flux
    Int viscosityModel;     // Flag for viscosity model. 0: Constant dynamic viscosity; 1: Sutherland's law

    Int SGSmodel;           // Flag for sub-grid scale (SGS) model. 0: No SGS model. 1: Smagorinsky/Yoshizawa/Knight model. 
                            //                                      2: WALE/Yoshizawa/Knight model. 3: Vreman/Yoshizawa/Knight model. SGS model only available for 3D solver.   
    Int convStabMethod;     // Flag for convective stabilization tensor. 0: Constant tau, 1: Lax-Friedrichs; 2: Roe.
    Int diffStabMethod;     // Flag for diffusive stabilization tensor. 0: No diffusive stabilization.
    Int rotatingFrame;      // Flag for rotating frame. Options: 0: Velocities are respect to a non-rotating frame. 1: Velocities are respect to a rotating frame.
};


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

// // element struct
// struct elementstruct {        
//     //elementinfo eleminfo;
//     Int dim;      // spatial dimension
//     Int elemtype; // type         
//     Int nfe; // number of faces    
//     Int nle;  // number of edges 
//     Int nve; // number of vertices 
//     vector< Int > facetype; // 0 triangular, 1 quadrilateral 
//     vector< vector< Int > > face; // face connectivity     
//     //vector< vector< Int > > edge; // edge connectivity     
//     
//     vector < Int > porder; // polynomial degrees along each direction
//     Int nodetype; // 0 uniform, 1 optimal        
//     Int npe; // number of nodes
//     Int nte; // number of subelements
//     vector< double>  plocvl; // nodal points
//     vector< Int > tlocvl;    // subelements
//     vector< Int > npf; // number of nodes on faces
//     vector< Int > ntf; // number of subelements on faces
//     vector< vector< double> > plocfc; // nodal points on faces
//     vector< vector< Int > > tlocfc;   // subelements on faces
//     vector< vector< Int > > perm;     // indices of the nodes on faces
// };

struct masterstruct {    
    Int dim;      // spatial dimension
    Int elemtype; // % element type  
    Int nodetype; // node type  
    Int npe; // number of DG nodes
    Int npv; // number of DG nodes
    Int nge; // number of Gauss nodes
    Int ngv; // number of Gauss nodes
    Int ngeR; // number of Gauss nodes
    Int ngvR; // number of Gauss nodes
    Int nfe; // number of faces    
    Int nle;  // number of edges 
    Int nve; // number of vertices 
    Int nplmax; // number number of DG nodes on each element vertice
    
    vector<Int> porder; // polynomial degrees along each direction
    vector<Int> pgauss;    
    vector<Int> pgaussR;    
    vector<Int> ndims;
    vector<Int> npf;
    vector<Int> ngf;
    vector<Int> ngfR;
    vector<Int> nvf;
    vector<Int> npl;
    
    vector<Int> permnode;
    vector<Int> permedge;
    vector< vector<Int> > perm;
    vector< vector<Int> > face;
    
    vector<double> plocvl;
    vector<Int> tlocvl;
    
    vector<double> gpvl;
    vector<double> gwvl;    
    vector<double> shapvl;
    vector<double> shapvt;
    vector<double> shapvg;
    vector<double> shapvgdotshapvl;
    vector<double> gpvlR;
    vector<double> gwvlR;
    vector<double> shapvtR;
    vector<double> shapvgR;
    
    vector< vector<double> > plocfc;
    vector< vector< Int > > tlocfc;  
    
    vector< vector<double> > gpfc;
    vector< vector<double> > gwfc;    
    vector< vector<double> > shapfc;
    vector< vector<double> > shapft;
    vector< vector<double> > shapfg;
    vector< vector<double> > shapfgdotshapfc;
    vector< vector<double> > gpfcR;
    vector< vector<double> > gwfcR;        
    vector< vector<double> > shapftR;
    vector< vector<double> > shapfgR;
    
    vector<double> shapnv;
    vector<double> shapnvt;   
    vector<double> projLowP;
        
    // NEW FIELDS FOR NEW MATRIX ASSEMBLY:
    //Int pgaussR;
    Int pgaussJ;
    Int pgaussQ;
    
    Int quadTypeR;
    Int quadTypeJ;
    Int quadTypeQ;
    
    Int nqvR;
    Int nqvQ;
    Int nqvJ;
    vector<double> gpvlJ;
    vector<double> gpvlQ;    
    vector<double> gwvlJ;
    vector<double> gwvlQ;    
    vector<double> shapvlJ;
    vector<double> shapvlQ;    
    vector<double> shapvtJ;
    vector<double> shapvtQ;    
    vector<double> shapvgJ;
    vector<double> shapvgQ;
    vector<double> shapvgdotshapvlR;
    vector<double> shapvgdotshapvlJ;
    vector<double> shapvgdotshapvlQ;
    
    vector<Int> nqfR;
    vector<Int> nqfQ;
    vector<Int> nqfJ;    
    vector< vector<double> > gpfcJ;
    vector< vector<double> > gpfcQ;    
    vector< vector<double> > gwfcJ;
    vector< vector<double> > gwfcQ;
    vector< vector<double> > shapfcJ;
    vector< vector<double> > shapfcQ;
    vector< vector<double> > shapftJ;
    vector< vector<double> > shapftQ;
    vector< vector<double> > shapfgJ;
    vector< vector<double> > shapfgQ;
    vector< vector<double> > shapfgdotshapfcR;
    vector< vector<double> > shapfgdotshapfcJ;
    vector< vector<double> > shapfgdotshapfcQ;
};
       
struct meshstruct {
    Int my_rank;
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
    
    vector<Int> nfes;
    vector<Int> npes;
    vector< vector<Int> > npfs;        
            
    vector<Int> ndims;            
    vector<double> p; // coordinates of vertice
    vector<double> dgnodes; // coordinates of DG nodes   
    vector<Int> t;  // element-to-vertice connectivities
    //vector<Int> dgnodesidx; // element index of dgnodes        
    vector<Int> elcon; // element-to-entity connectivities
    //vector<Int> elconidx; // element index of elcon        
    //vector< vector<Int> > npfidx;  // face index of npf
    vector<Int> bf;  // numbering of boundary faces
    vector<Int> t2f; // element-to-face connectivities
    vector<Int> t2t; // element-to-element connectivities
    vector<Int> f;   // face-to-element connnectivities 
    vector<Int> extintelem;
    vector<Int> extintface;
    vector<Int> extintent;
    vector<Int> elem2cpu;
    vector<Int> face2cpu;
    vector<Int> ent2cpu;
    vector<Int> nbsd;
    vector<Int> cg2dg;
    vector<Int> dg2cg;             
            
    vector<Int> porder;    
    vector<Int> pgauss;
    vector<Int> quadtype;    
    //vector<elementstruct> elements; // types of elements in the mesh    
    //vector<masterstruct> masters; // types of masters in the mesh    
    vector<Int> elementtype;  // associated with the structure elements 
    vector<Int> physics;     // associated with the governing equations            
    vector<Int> isEDGface;    
    
    vector<double> elemMeasure; // Element measure: ne
    vector<double> hAvg;        // Characteristic element size: ne
    vector<double> M;           // Metric tensor at dgnodes
    vector<double> Minv;        // Inverse of metric tensor at dgnodes
};

struct dmdstruct {    
    Int my_rank; // index of the processor
    vector<Int> nbsd;    // neighboring cpus
    vector<Int> intelem; // nonoverlapping global elements
    vector<Int> intent;  // nonoverlapping global entities
    vector<Int> elempart; // overlapping global elements
    vector<Int> elempartpts; // classifiers of overlapping global elements: (interior, interface, exterior)
    vector<Int> entpart;  // overlapping global entities
    vector<Int> entpartpts; // classifiers of overlapping global entities: (interior, interface, exterior)
    vector<Int> elemrecv;  // local elements received from neighboring cpus
    vector<Int> elemrecvpts; // classifiers of local elements received from neighboring cpus: (# elements from ncpu1, ...) 
    vector<Int> elemsend;  // local elements sent to neighboring cpus
    vector<Int> elemsendpts; // classifiers of local elements sent to neighboring cpus: (# elements to ncpu1, ...)    
    vector<Int> entrecv;    // local entities received from neighboring cpus
    vector<Int> entrecvpts; // classifiers of local entities received from neighboring cpus: (# entities from ncpu1, ...)          
    vector<Int> entsend;  // local entities sent to neighboring cpus
    vector<Int> entsendpts; // classifiers of local entities sent to neighboring cpus: (# elements to ncpu1, ...)    
    vector<Int> vecrecv;  // local vectors received from neighboring cpus to perform matrix-vector product
    vector<Int> vecrecvpts; // classifiers of local vectors received from neighboring cpus: (# vectors from ncpu1, ...)          
    vector<Int> vecsend;    // local vectors sent to neighboring cpus to perform matrix-vector product
    vector<Int> vecsendpts; // classifiers of local vectors sent to neighboring cpus: (# vectors to ncpu1, ...)             
    vector<Int> matrecv;    // local matrices received from neighboring cpus to construct the preconditioner
    vector<Int> matrecvpts; // classifiers of local matrices received from neighboring cpus: (# matrices from ncpu1, ...)                
    vector<Int> matsend;   // local matrices sent to neighboring cpus to construct the preconditioner
    vector<Int> matsendpts;  // classifiers of local matrices sent to neighboring cpus: (# matrices to ncpu1, ...)                                         
    vector<Int> rowent2elem; // global entity-to-element connectivities
    vector<Int> colent2elem; // global entity-to-element connectivities
    vector<Int> rowent2ent;  // global entity-to-entity connectivities
    vector<Int> colent2ent;  // global entity-to-entity connectivities
    vector<Int> bcrs_rowent2elem; // local entity-to-element connectivities
    vector<Int> bcrs_colent2elem; // local entity-to-element connectivities
    vector<Int> bcrs_rowent2ent;  // local entity-to-entity connectivities
    vector<Int> bcrs_colent2ent;  // local entity-to-entity connectivities    
    vector<Int> ent2ind;  // global-to-local entity mapping  
    vector<Int> elcon;    // local element-to-entity connectivities
    vector<Int> t2f;      // local elemeent-to-face connectivities
    vector<Int> t2t;      // local elemeent-to-element  connectivities
    Int  maxBlocksPerRow; // maximum number of entities per row
    Int  minBlocksPerRow; // minimum number of entities per row    
};

struct solstruct {
    vector<Int> ndims;
    vector<double> UDG; /* UDG = (U, Q, P) */
    vector<double> UH;  /* UH = UHAT */
    vector<double> PDG; /* additional unknown for wave applications */
    vector<double> SH;  /* Source from the discretization of the time derivatives */
    vector<double> SP;  /* additional source for wave applications */

    vector<double> UDG4FD; /* Only used for matrix-free matrix-vector product */
    vector<double> UH4FD;  /* Only used for matrix-free matrix-vector product */
    
    vector<double> UDG_avg; // Time average UDG
    vector<double> UH_avg;  // Time average UH
    
    vector<double> Un; /* Un = (U, UHAT) */
    vector<double> Um; /* Um = (U, UHAT) */
    vector<double> Vn; /* Vn = (U, UHAT) */
    vector<double> R;  /* R = (RU, RUHAT) */

    vector<double> avField_DG;  // DG artificial viscosity field
    vector<double> avField_p1CG;  // p=1 CG artificial viscosity field
    vector<double> avField;     // NEW FIELD FOR NEW MATRIX ASSEMBLY

    vector<double> UDG2Write; // Auxiliary variable to write solution to a binary file. Only required if writeQflag = 1.

    // For Minimal Residual algorithm to compute the initial guess for Newton iteration
    vector<double> UDG_initMR;
    vector<double> UH_initMR;
    
    /* for time-dependent problems only */
    vector<double> UDGi; /* store UDG from the previous timestep */
    vector<double> UHi;  /* store UH from the previous timestep */
    vector<double> PDGi; /* store PDG from the previous timestep */
    
    /* used to recover DU */
    vector<double> DinvRu;  /* inverse(D)*Ru */
    vector<double> DinvF;   /* inverse(D)*F */

    /* Only needed for Quasi Newton. TODO: Wouldn't it be better to store F instead of DinvF for quasi-Newton? */
    vector<double> Dinv;    /* inv( dRu_bar/du ) */
    vector<float> DinvFloat;    /* inv( dRu_bar/du ) */
    vector<Int> ipivD;      /* pivots for LU factors of dRu_bar/du */
    vector<double> K;       /* dRh/du */
    vector<float> Kfloat;       /* dRh/du */

    double alpha;  /* stepsize for damped Newton */
    double rNorm;  /* the norm of the residual vector for both U and UH */
    
    vector<double> Cs;      // Constant for dynamic Smagorinsky SGS model: npv / ne
    
    /* call to write (UDG,UH) into a binary file */
    void writeSol2File(string filename, Int ne, Int n, Int* ndims, Int writeQflag);
    
    /* call to write (UDG_avg,UH_avg) into a binary file */
    void writeAvgSol2File(string filename_avg, Int ne, Int n, Int* ndims, Int writeQflag);
};

//void mkmesh(meshstruct &mesh, vector<Int> &porder, Int dim, Int nodetype);

#endif
