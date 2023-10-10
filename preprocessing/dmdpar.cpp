#ifndef __DMDPAR
#define __DMDPAR

// Written by: C. Nguyen & P. Fernandez

#include "mkdmd.cpp"

//dmdparallel(sys, mesh, app, sol, dmd);            
void dmdparallel(sysstruct &sys, meshstruct &mesh, appstruct &app, solstruct &sol, dmdstruct &dmd)
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
        if (sys.my_rank==0)
            printf("HDG domain decomposition...\n");     
        domaindecomposition(dmd, mesh.t2f, mesh.elcon, mesh.t2f, mesh.elementtype, mesh.extintelem, mesh.extintface,
                mesh.elem2cpu, mesh.face2cpu, mesh.nfes, mesh.npfs, nie, nin, nxe, nxn, 1, mesh.my_rank);                                
        MPI_Barrier(MPI_COMM_WORLD);
        Int nnbsd = (Int) dmd.nbsd.size();                  
        vector< dmdstruct > nbdmd(nnbsd, dmdstruct());
        if (sys.my_rank==0)
           printf("Do make_nbdmd...\n");                       
        make_nbdmd(nbdmd, dmd);     
        MPI_Barrier(MPI_COMM_WORLD);
        if (sys.my_rank==0)
            printf("Do sendrecvhdg...\n");     
        sendrecvhdg(dmd);
        MPI_Barrier(MPI_COMM_WORLD);
        if (sys.my_rank==0)
            cout<<sys.my_rank<<" Finish sendrecvhdg"<<endl;
    }
    else {                  
        if (sys.my_rank==0)
            cout<<"EDG domain decomposition..."<<endl;     
        domaindecomposition(dmd, mesh.elcon, mesh.elcon, mesh.t2f, mesh.elementtype, mesh.extintelem, mesh.extintent,
               mesh.elem2cpu, mesh.ent2cpu, mesh.nfes, mesh.npfs, nie, nin, nxe, nxn, 0, mesh.my_rank);                                
        //error("Finish domain decomposition\n");
        MPI_Barrier(MPI_COMM_WORLD);  
        if (sys.my_rank==0)
            cout<<"Finish EDG domain decomposition"<<endl;
        Int nnbsd = (Int) dmd.nbsd.size();   
        vector< dmdstruct > nbdmd(nnbsd, dmdstruct());
        if (sys.my_rank==0)
            cout<<"Do make_nbdmd..."<<endl;                  
        make_nbdmd(nbdmd, dmd);                
        MPI_Barrier(MPI_COMM_WORLD);
        if (sys.my_rank==0)
            cout<<"Finish mk_nbdmd"<<endl;
        if (sys.my_rank==0)
            cout<<"Do sendrecvedg..."<<endl;     
        sendrecvedg(dmd, nbdmd);
        MPI_Barrier(MPI_COMM_WORLD);  
        if (sys.my_rank==0)
            cout<<"Finish sendrecvedg"<<endl;
    }                        
    
    // update mesh
    mesh.ne = dmd.elemmap.size();    
    mesh.nf = *max_element(dmd.t2f.begin(), dmd.t2f.end())+1;
    mesh.ndofuh = *max_element(dmd.elcon.begin(), dmd.elcon.end())+1;
    
    mesh.t2f = dmd.t2f; //iarray2datn(dmd.t2f, mesh.ne, mesh.nfemax);
    mesh.elcon = dmd.elcon; //iarray2datn(dmd.elcon, ne, mesh.nfemax*mesh.npfmax);
    
    vector<double> tm;
    vector<Int> tn;
    tm = darray2datindex(sol.UDG, dmd.elemmap, mesh.npemax*app.nc); sol.UDG = tm;
    tm = darray2datindex(mesh.dgnodes, dmd.elemmap, mesh.nmemax*app.ncd); mesh.dgnodes = tm;           
    tn = iarray2datindex(mesh.bf, dmd.elemmap, mesh.nfemax); mesh.bf = tn;    
    tn = iarrayatindex(mesh.elementtype, dmd.elemmap); mesh.elementtype = tn;
    tn = iarrayatindex(mesh.isEDGelement, dmd.elemmap); mesh.isEDGelement = tn;
    tn = iarray2datindex(mesh.t, dmd.elemmap, mesh.nvemax); mesh.t = tn;  
    tn = iarray2datindex(mesh.t2t, dmd.elemmap, mesh.nfemax); mesh.t2t = tn;   
    tn = iarrayatindex(mesh.elem2cpu, dmd.elemmap); mesh.elem2cpu = tn;
    tn = iarrayatindex(mesh.extintelem, dmd.elemmap); mesh.extintelem = tn;   
            
    if (hybrid == 0)  { // HDG 
        tm = darray2datindex(sol.UH, dmd.entmap, app.nch*mesh.npfmax); 
        sol.UH = tm;   
    }
    else {
        tm = darray2datindex(sol.UH, dmd.entmap, app.nch); 
        sol.UH = tm;         
    }
    
    // update sys    
    sys.entpart = dmd.entpart;
    sys.entpartpts = dmd.entpartpts;
    sys.entrecv = dmd.entrecv;
    sys.entrecvpts = dmd.entrecvpts;
    sys.entsend = dmd.entsend;
    sys.entsendpts = dmd.entsendpts;
    sys.elempart = dmd.elempart;
    sys.elempartpts = dmd.elempartpts;
    sys.elemrecv = dmd.elemrecv;
    sys.elemrecvpts = dmd.elemrecvpts;
    sys.elemsend = dmd.elemsend;
    sys.elemsendpts = dmd.elemsendpts;
    sys.vecrecv = dmd.vecrecv;
    sys.vecrecvpts = dmd.vecrecvpts;
    sys.vecsend = dmd.vecsend;
    sys.vecsendpts = dmd.vecsendpts;    
    sys.matrecv = dmd.matrecv;
    sys.matrecvpts = dmd.matrecvpts;
    sys.matsend = dmd.matsend;
    sys.matsendpts = dmd.matsendpts;    
    sys.ent2entStart = dmd.bcrs_rowent2ent;
    sys.ent2ent = dmd.bcrs_colent2ent;
    sys.nbsd = dmd.nbsd;
    
    sys.nproc = app.nproc;
    sys.my_rank = app.my_rank;    
    if (hybrid == 0) //HDG
        sys.blkSize = app.nch*mesh.npfmax;
    else
        sys.blkSize = app.nch;    
    sys.nentrecv = sys.entrecv.size();
    sys.nentsend = sys.entsend.size();    
    sys.nelemrecv = sys.elemrecv.size();
    sys.nelemsend = sys.elemsend.size();
    sys.nvecrecv = sys.matrecv.size();
    sys.nvecsend = sys.matsend.size();
    sys.nmatrecv = sys.matrecv.size();
    sys.nmatsend = sys.matsend.size();    
    sys.nnbsd = sys.nbsd.size();    
    sys.BJ_nrows = sys.entpartpts[0]+sys.entpartpts[1];
    sys.numEntities = sys.entpart.size();
    sys.numBlocks = sys.ent2ent.size();    
    
    Int a;
    sys.maxBlocksPerRow = 0;
    sys.minBlocksPerRow = 9999999;
    for (Int i=1; i<sys.ent2entStart.size(); i++) {
        a = sys.ent2entStart[i]-sys.ent2entStart[i-1];
        if (a > sys.maxBlocksPerRow)
            sys.maxBlocksPerRow = a;
        if (a < sys.minBlocksPerRow)
            sys.minBlocksPerRow = a;
    }
}    

#endif
