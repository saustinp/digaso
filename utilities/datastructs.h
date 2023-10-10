#ifndef __DATASTRUCTS_H
#define __DATASTRUCTS_H

// Written by: C. Nguyen & P. Fernandez

using namespace std;

struct solstruct {
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

struct elemstruct {
/* storing the following system for one element
 [M -C  E] [dq] = [Rq]
 [B  D  F] [du] = [Ru]
 [G  K  H] [dh] = [Rh]
*/
    double*  M;
    double*  E;
    double*  C;
    double*  Rq;
    double*  BD;
    double*  F;
    double*  Ru;
    double*  GK;
    double*  H;
    double*  Rh;
    double*  Rhonly;
    vector<double> D1;
    vector<double> F1;
    vector<double> K1;
    vector<double> H1;
    vector<double> K;

    Int my_rank;
    
    vector<double> F_dense;             ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> F_conv;              ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    double*  D_tmp;                     ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> K_tmp;               ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    double*  D_LU;                      ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    double*  D_LU_tmp;                  ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> D_inv;               ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> D_inv_extCol;        ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> DiF_conv;            ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> DiF_tmp;             ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> DiF_extRow;          ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    double*  DiF_ij;                    ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> Di_tmp;              ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    double*  Ru_i;                      ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<Int> pivDii;                 ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> workDii;             ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> workD;               ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
};

struct sysstruct {
    /* TODO: Create an structure for the preconditioner */
    Int numEntities;        // Number of entities in the mesh (nfn in EDG or nf in HDG)
    Int numBlocks;          // Number of blocks in matrix of linear system
    Int blkSize;            // Block size (nch in EDG or nch*npf in HDG)
    Int maxBlocksPerRow;    // Maximum number of blocks per row in the matrix of the linear system
    Int BJ_nrows;           // Number of rows associated to entities in the processor
    Int BJ_nblks;           // Number of blocks in rows associated to entities in the processor
    Int BK_nrows;           // Number of rows associated to entities in the processor that contain columns not in the processor
    Int BK_nblks;           // Number of blocks in rows associated to entities in the processor that contain columns not in the processor
    Int nproc;              // Number of processors
    Int nentpartpts;        
    Int nentrecv;
    Int nentsend;
    Int nelempartpts;        
    Int nelemrecv;
    Int nelemsend;
    Int nnbsd;
    Int my_rank;            // Index of the processor
    Int noThreads;         // Number of OpenMP threads
    Int computeGlobalEnt2entWeight;       // Only applies for MPI code. If set to 1, the Frobenius norm of the Schur matrix blocks are computed, stored in a file and the execution is terminated.
    Int nmatrecv;
    Int nmatsend;
    
    vector<Int> ent2ent;
    vector<Int> ent2entStart;
//    vector<Int> globalEnt2ent;      // Only used in (1) MPI code and (2) if sys.computeGlobalEnt2entWeight == 1
//    vector<Int> globalEnt2entStart; // Only used in (1) MPI code and (2) if sys.computeGlobalEnt2entWeight == 1
    vector<Int> oent2ent;           // Only required for implementation No. 2 of BILU0
    vector<Int> oent2entStart;      // Only required for implementation No. 2 of BILU0
    vector<Int> LUoent2ent;         // Only required for implementation No. 3 of BILU0
    vector<Int> Loent2entStart;     // Only required for implementation No. 3 of BILU0
    vector<Int> Uoent2entStart;     // Only required for implementation No. 3 of BILU0
    vector<Int> BJ_rowpts;
    vector<Int> BJ_colind;
    vector<Int> BK_rowpts;
    vector<Int> BK_colind;
    vector<Int> BK_rowind;
    vector<double> ent2entWeight;       // Only from entities allocated to the processor (weights from entities in other processors are treated by those other processors)
//     vector<Int> ent2entWeightLen;       // Length of the vector to be received by processor No. 1 from each processor
    vector<Int> entpart;                // local-to-global entity mapping
    vector<Int> entpartpts;
    vector<Int> entrecv;
    vector<Int> entrecvpts;
    vector<Int> entsend;
    vector<Int> entsendpts;
    vector<Int> elempart;
    vector<Int> elempartpts;
    vector<Int> elemrecv;
    vector<Int> elemrecvpts;
    vector<Int> elemsend;
    vector<Int> elemsendpts;
    vector<Int> nbsd;
    vector<Int> matrecv;
    vector<Int> matrecvpts;
    vector<Int> matsend;
    vector<Int> matsendpts;
    
    // Variables for Minimal Residual algorithm for non-linear initial guess:
    vector<double> Ru_MR;
    vector<double> dRuda_MR;
    vector<double> a2_MR;
    vector<vector<double> > a2s_MR;
    vector<double> da_MR;
    vector<double> C_MR;
    vector<double> Clocal_MR;
    
    Int orthogMethod;               // Orthogonalization method in GMRES. 0: MGS, 1: ICGS, 2: IMGS (only MGS and ICGS are available in MPI code), 3: ICGS without convergence guarantee
    Int maxiter;                    // Maximum number of GMRES iterations
    Int restart;                    // Parameter k in GMRES(k)
    double tol;                     // Convergence tolerance for linear system
    Int preconditionerSide;         // 0: Left preconditioner. 1: Right preconditioner
    Int reorderMethod;              // 0: No reordering. 1: Approximate (inexact) MDF. 2: Exact MDF. 3: Approximate (inexact) MDF with constraints (only valid for parallel BJ preconditioner). 4: Exact MDF with constraints (only valid for parallel BJ preconditioner)
    Int preconditioner;             // -1: No preconditioner, 0: Restricted Additive Schwarz (RAS) + BILU(k), 1: Subdomain-wise Block Jacobi (BJ) + BILU(0), 2: Entity-wise block Jacobi
    Int NewtonMaxiter;              // Maximum number of Newton/quasi-Newton iterations
    Int trueNewtonMaxiter;          // Maximum number of true Newton iterations
    double NewtonTol;               // Newton tolerance for non-linear system

    double rNorm;                   /* the norm of the residual vector for both U and UH */

    Int print;                      // 0: no print, 1: print on screen
    Int schurImplementation;        // Implementation flag. 0 and 1 values are valid
    Int matvecImplementation;       // Implementation flag for matrix-vector product. Options: 0, 1, 2 (finite differences for Gateaux derivative)
    Int precSolveImplementation;    // Implementation flag for preconditioner solve. Options: 0 and 1
    Int precPrecision;              // Precision for preconditioner solve. Options: 0: Single precision. 1: Double precision
    Int matvecPrecision;            // Precision for matrix-vector product. Options: 0: Single precision. 1: Double precision
    Int orthogPrecision;            // Precision for orthogonalization. Options: 0: Single precision. 1: Double precision
    Int quasiNewtonAccuracy;        // Format in which Dinv and K are stored for quasi-Newton. 0: Single precision; 1: Double precision
    Int adaptiveGMREStol;           // Adaptive GMRES tolerance for nonlinear problems. Options: 1: Adaptive method. 0: Non-adaptive method. Note the adaptive strategy is not used for linear PDEs regardless the value of this flag.
    Int adaptGMREStol;              // Only applies for linear problems with adaptiveGMREStol == 1. Has two meanings:
                                    // Before executing GMRES routine: If adapting GMRES tolerance is allowed.
                                    // After executing GMRES routine: If GMRES tolerance was adapted in the previous linear solve
    Int robustMode;                 // Options: 0 (default): All flags to default value. 1: All flags to safety (or robust) mode
    
    Int linearProblem;      /* flag to determine linear or nonlinear problems */
    long linearSolvesClocks;         // Clocks spent in lear solves since the Jacobian matrix was computed the last time (only applies for quasi-Newton)
    long lastAssemblyAndPrecClocks;  // Clocks spent to perform the last matrix assembly and preconditioner computation (only applies for quasi-Newton)

    vector<double> Hg;              // Part of global matrix such that M_i ~= (Hg_i)^-1 (BILU0 inverse)
    vector<double> Kg;              // Part of global matrix not contained in Hg
    vector<double> Rg;              // Global RHS (this vector will be modified during GMRES)
    vector<double> Rg4FD;           // Global RHS (used for matrix-free matrix-vector product)
    vector<double> Rg_0;            // Global RHS (this vector will not be modified during GMRES)
    vector<double> x;               // Solution to linear system (Hg * x = Rg)
    vector<vector<double> > xDenseRow;       // Dense version of x for a particular row
//    vector<double> xDense;          // Dense version of x
    vector<double> Mg;              // Preconditioner matrix
    vector<double> r;               // Residual vector
    vector<double> r4FD;            // Residual vector (used for matrix-free matrix-vector product)
    vector<double> v;               // Krylov vectors (also have other uses)

    vector<double> s;               // Auxiliary variable for orthogonalization
    vector<double> stmp;            // Auxiliary variable for orthogonalization
    vector<vector<double> > stmps;            // Auxiliary variable for orthogonalization
    vector<double> Mx;              // Preconditioner matrix apply to the solution vector
    double* Mv;                     // Preconditioner matrix apply to the last vector in the Arnoldi basis

    // Single precision arrays for mixed-precision algorithms:
    vector<float> Hg_sp;
    vector<float> Kg_sp;
    vector<float> Mg_sp;
    vector<float> x_sp;
    vector<vector<float> > xDenseRow_sp;       // Dense version of x_sp for a particular row
    vector<float> v_sp;            // Krylov vectors in single precision
    vector<float> s_sp;               // Auxiliary variable for orthogonalization
    vector<float> stmp_sp;            // Auxiliary variable for orthogonalization
    vector<vector<float> > stmps_sp;            // Auxiliary variable for orthogonalization
    vector<float> r_sp;
    vector<float> rtmp_sp;
    vector<vector<float> > rtmps_sp;
    vector<float> rDenseRow_sp;

    vector<Int> ipiv;                   // Pivots for LU factors of diagonal blocks of BILU(k) preconditioner
    vector<vector<Int> > ipivs;
    vector<Int> ordered2unordered;      // Mapping from ordered entities to unordered entities
    vector<Int> unordered2ordered;      // Mapping from unordered entities to ordered entities

    double* C;
    double* w;

    vector<double> e1;
    vector<double> rev;
    vector<double> y;
    vector<double> hy;
    vector<double> hh;
    vector<double> hhls;
    vector<double> rtmp;
    vector<vector<double> > rtmps;
    vector<double> work;
    vector<vector<double> > works;
    vector<double> workls;

    vector<double> buffsend;
    vector<double> buffrecv;
    vector<double> buffrecvEnt2entWeight;
    vector<double> buffsendmat;
    vector<double> buffrecvmat;
    
    // Auxiliary vectors for approximate MDF ordering:
    vector<double> HiiInv;
    vector<vector<double> > HiiInvs;
    vector<double> HiiInvHij;
    vector<vector<double> > HiiInvHijs;
    
    // Auxiliary vectors for exact MDF ordering:
    vector<double> HkkInv;
    vector<vector<double> > HkkInvs;
    vector<double> HikHkkInv;
    vector<vector<double> > HikHkkInvs;
    vector<double> HikHkkInvHkj;
    vector<vector<double> > HikHkkInvHkjs;
    
#ifdef  HAVE_MPI
    MPI_Request * requests;
    MPI_Status * statuses;

//    MPI_Request * requestsEnt2entWeight;
//    MPI_Status * statusesEnt2entWeight;
#endif   

};


struct meshstruct {
    vector<double> dgnodes;
    vector<Int> elcon;
    vector<Int> elconintf;
    vector<Int> bf;
    vector<Int> t2f;
    vector<Int> t;
    vector<Int> perm;
    vector<Int> permgeom;
    vector<Int> cg2dg;
    vector<Int> dg2cg;
    vector<Int> isEDGelement;
    Int elemtype;
    
    vector<double> elemMeasure; // Element measure: ne
    vector<double> hAvg;        // Characteristic element size: ne
    vector<double> M;           // Metric tensor at dgnodes
    vector<double> Minv;        // Inverse of metric tensor at dgnodes
};


struct masterstruct {
    Int nd;
    Int porder;
    Int pgauss;
    Int npv;
    Int npf;
    Int ngv;
    Int ngf;
    Int elemtype;
    Int nodetype;
    
    vector<double> plocvl;
    vector<double> gpvl;
    vector<double> gwvl;
    vector<double> plocfc;
    vector<double> gpfc;
    vector<double> gwfc;
    vector<double> shapvl;
    vector<double> shapvt;
    vector<double> shapvg;
    vector<double> shapvgdotshapvl;
    vector<double> shapfc;
    vector<double> shapft;
    vector<double> shapfg;
    vector<double> shapfgdotshapfc;
    vector<double> shapmv;
    vector<double> shapmf;
    vector<Int> perm;
    
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
    
    Int nqfR;
    Int nqfQ;
    Int nqfJ;
    vector<double> gpfcR;
    vector<double> gpfcJ;
    vector<double> gpfcQ;
    vector<double> gwfcR;
    vector<double> gwfcJ;
    vector<double> gwfcQ;
    vector<double> shapfcR;
    vector<double> shapfcJ;
    vector<double> shapfcQ;
    vector<double> shapftR;
    vector<double> shapftJ;
    vector<double> shapftQ;
    vector<double> shapfgR;
    vector<double> shapfgJ;
    vector<double> shapfgQ;
    vector<double> shapfgdotshapfcR;
    vector<double> shapfgdotshapfcJ;
    vector<double> shapfgdotshapfcQ;
};


struct appstruct {
    vector<Int> bcm;
    vector<double> bcs;
    vector<Int> bcd;
    vector<double> bcv;
    vector<double> dt;      // TODO: Allow scalar field for local time stepping in steady problems?
    vector<double> param;
    vector<Int> flag;
    vector<double> factor;
    vector<Int> problem;
    Int porder;             // Polynomial order of numerical approximation
    Int my_rank;
    
    string filein;          // Name of binary file with input data
    string fileout;         // Name of binary file to write the solution

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

    vector<Int> DIRKnotConverged;
    vector<Int> BDFnotConverged;
    
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

struct tempstruct {
    double* pg;
    double* Jg;
    double* jacg;
    double* Xxg;

    double* pf;
    double* pgf;
    double* Jgf;
    double* nlgf;
    double* jacgf;

    double* nlgjac;
    double* udgg;
    double* udgg_ref;
    double* f;
    double* f_udg;
    double* s;
    double* s_udg;
    double* wrk;
    double* wrl;

    double* uh;
    double* uh_ref;
    double* uf;
    double* uf_ref;
    double* ugf;
    double* ugf_ref;
    double* uhg;
    double* uhrefg;
    double* avg_p1CG;
    double* avf_p1CG;
    double* avfg_p1CG;
    double* fh;
    double* fh_u;
    double* fh_uh;
    double* Rutmp;
    double* BDtmp;
    double* BDt;    
    double* Ft;
    double* Ftmp;
    double* Etmp;
    double* fb;
    double* fb_u;
    double* fb_uh;
    double* pft;
    double* uft;
    double* uht;
    double* nlt;
    double* uft_ref;
    double* uht_ref;
    double* avft_p1CG;

    double* BMiE;
    double* BMiC;
    double* GMiE;
    double* GMiC;        
    
    vector<double> shapMiC; //(ngv*npv*nd); 
    vector<double> shapMiE; //(ngv*ndf*nd); 
    vector<double> wrlshapMiC; //(ngv*nd1*ncu*npv*ncu);        
    vector<double> wrlshapMiE; //(ngv*nd1*ncu*ncu*ndf);
    
    vector<double> MiCf; 
    vector<double> shapMiCf; 
    vector<double> fhushapMiCf;
    vector<double> Df;
    
    vector<double> MiEf; 
    vector<double> shapMiEf; 
    vector<double> fhushapMiEf;
    vector<double> Ff;
    
    vector<Int> ipiv;
    vector<Int> ind;
    vector<Int> jnd;

    vector<double> H_tmp;
    vector<double> K_tmp;
    vector<double> F_tmp;
    vector<double> Rh_tmp;
    
    
    
    
    
    
    
    
    // NEW FIELDS FOR NEW MATRIX ASSEMBLY:
    double* pp;
    double* udgp;
    double* udgp_ref;
    double* avp;
    
    double* avf;
    
    double* fhR;
    double* fhR_u;
    double* fhR_uh;
    double* uhR;
    double* ufR;
    double* pfR;
    double* nlR;
    double* jacfR;
    double* nljacR;//NOT NECESSARY
    double* uh_refR;
    double* uf_refR;
    double* avfR;
    
    double* fhJ;
    double* fhJ_u;
    double* fhJ_uh;
    double* uhJ;
    double* ufJ;
    double* pfJ;
    double* nlJ;
    double* jacfJ;
    double* nljacJ;
    double* uh_refJ;
    double* uf_refJ;
    double* avfJ;
    
    double* fhjacR;
    double* fh_ujacJ;
    double* fh_uhjacJ;
    
    double* fbR;
    double* fbR_u;
    double* fbR_uh;
    double* uftR;
    double* uhtR;
    double* uft_refR;
    double* uht_refR;
    double* pftR;
    double* JfR;
    double* XxfR;
    double* nltR;
    double* avftR;
    
    double* fbJ;
    double* fbJ_u;
    double* fbJ_uh;
    double* uftJ;
    double* uhtJ;
    double* uft_refJ;
    double* uht_refJ;
    double* pftJ;
    double* JfJ;
    double* XxfJ;
    double* nltJ;
    double* avftJ;
    
    double* shapMiCfJ;
    double* fh_ushapMiCfJ;
    double* shapMiEfJ;
    double* fh_ushapMiEfJ;
    
    double* pvQ;
    double* Jv_Q;
    double* jacvQ;
    double* XxvQ;
    
    double* pfQ;
    double* Jf_Q;
    double* XxfQ;
    double* jacfQ;
    double* nlfQ;
    double* nlfjacQ;
    
    double* sJ;
    double* fJ;
    double* sJ_udg;
    double* fJ_udg;
    double* jacJ;
    double* XxJ;
    double* pJ;
    double* J_J;
    
    double* sR;
    double* fR;
    double* sR_udg;      // NOT NECESSARY
    double* fR_udg;      // NOT NECESSARY
    double* jacR;
    double* XxR;
    double* pR;
    double* J_R;
    
    double* shapMiC_J;
    double* shapMiE_J;
    double* wrlshapMiC_J;
    double* wrlshapMiE_J;
};

// struct setupstruct {
//     /* Non-linear solver setup */
//     int maxNRiter;
//     double NRtol;
//     double initialAlpha;
//     double contractionAlpha;
//     double minAlpha;
//
//     /* GMRES setup */
//     int restart;
//     int maxRestarts;
//     int maxGMRESiter;
//     double GMREStol;
//     int orthogMethod;
//
//     /* Preconditioner setuip */
//     int preconditionerSide; /* 0: left preconditioner; 1: right preconditioner */
// };

#endif

