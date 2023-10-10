#ifndef __PREPROCESSING
#define __PREPROCESSING

#include "errormsg.cpp"
#include "mathutilities.cpp"
#include "../utilities/block_crs_k.cpp"
#include "../gmres/computeMgMinusHg.cpp"
#include "ioutilities.cpp"
#include "readbinaryfiles.cpp"
#include "initializestructs.cpp"
#include "../master/mkshape.cpp"
//#include "mkdmd.cpp"

#ifdef HAVE_MPI
void dmdparallel(sysstruct &sys, meshstruct &mesh, appstruct &app, solstruct &sol, dmdstruct &dmd);
#endif

#ifndef HAVE_MPI              
void dmdserial(sysstruct &sys, meshstruct &mesh, appstruct &app);
#endif

void preprocessing(meshstruct &mesh, vector< masterstruct > &master, appstruct &app, solstruct &sol, elemstruct &elems, sysstruct &sys, tempstruct &temps, string filein, string fileout) 
{
    // read from binary files to construct app, master, mesh, sol
    if (sys.my_rank==0)
        printf("Read data files...\n");         
    readInput(app, master, mesh, sol, filein);                    
    if (sys.my_rank==0)
        printf("End reading data files...\n");         
    
    if (app.nproc>1) {
#ifdef HAVE_MPI       
        if (sys.my_rank==0)
            printf("Create domain decomposition...\n");     
        dmdstruct dmd;
        dmdparallel(sys, mesh, app, sol, dmd);            
        MPI_Barrier(MPI_COMM_WORLD);
        if (sys.my_rank==0)
            cout<<"Finish domain decomposition"<<endl;
            
//         if (sys.my_rank==0)
//             printf("Write solution information file...");     
//                 
//         if (app.debugmode == 1) 
//             writeOutput(app, master, mesh, sol, dmd, fileout);
// 
//             //vector<Int> ndims = getndims(mesh, master, app, sol, sys);
//                //printiarray(ndims);             
//         Int N = dmd.elempartpts[0]+dmd.elempartpts[1];
//         string filename1 = fileout + "elempart" + ".bin";
//         mpiwriteiarray2file(filename1, &dmd.elempart[0], N, app.my_rank, app.nproc);                
// 
//         N = dmd.entpartpts[0]+dmd.entpartpts[1];
//         string filename2 = fileout + "entpart" + ".bin";
//         mpiwriteiarray2file(filename2, &dmd.entpart[0], N, app.my_rank, app.nproc);

        Int N = dmd.elempartpts[0]+dmd.elempartpts[1];
        string filename1 = fileout + "elempart" +  "_np" + NumberToString(app.my_rank) + ".bin";            
        writeiarray2file(filename1, &dmd.elempart[0], N);
        
        N = dmd.entpartpts[0]+dmd.entpartpts[1];
        string filename2 = fileout + "entpart" +  "_np" + NumberToString(app.my_rank) + ".bin";            
        writeiarray2file(filename2, &dmd.entpart[0], N);
        
        // clear memory
        cleardmdstruct(dmd); 
        
        MPI_Barrier(MPI_COMM_WORLD);
        if (sys.my_rank==0)
            cout<<"End preprocessing"<<endl;                
#endif                
    }
    else {
#ifndef HAVE_MPI              
        dmdserial(sys, mesh, app);    
#endif    
    }
    
    // allocate memory     
    if (sys.my_rank==0)
        printf("Allocate memory...\n");     
    initializeStructs(mesh, master, app, sol, elems, sys, temps);   
    
    // Compute additional fields in master structure (required for the metric M and the inverse metric Minv tensor fields)
    for (int iM = 0; iM < master.size(); iM++) {
        if (master[iM].npe > 0) {
            Int npe = master[iM].npe;
            Int nd = master[iM].nd;
            Int porder = master[iM].porder[0];
            Int elemtype = master[iM].elemtype;
            master[iM].shapnv.resize(npe*npe*(nd+1));
            mkshape(&master[iM].shapnv, &master[iM].plocvl[0], &master[iM].plocvl[0], npe, elemtype, porder, nd, npe);
            master[iM].shapnvt.resize(npe*npe*(nd+1));
            for (int i = 0; i < (nd+1); i++)
                for (int j = 0; j < npe; j++)
                    for (int k = 0; k < npe; k++)
                        master[iM].shapnvt[i*npe*npe+k*npe+j] = master[iM].shapnv[i*npe*npe+j*npe+k];
        }
    }
}

// void dmdserial(sysstruct &sys, meshstruct &mesh, appstruct &app)
// {
//     Int N = mesh.ndims.size();
//     Int hybrid = mesh.ndims[N-2];        
//     
//     mesh.ne = mesh.elementtype.size();    
//     mesh.nf = *max_element(mesh.t2f.begin(), mesh.t2f.end())+1;
//     mesh.ndofuh = *max_element(mesh.elcon.begin(), mesh.elcon.end())+1;
//     
//     vector<Int> rowent2elem, colent2elem;
//     if (hybrid == 0) { // HDG                
//         elcon2entcon(rowent2elem,colent2elem,sys.ent2entStart,sys.ent2ent,mesh.t2f,mesh.elementtype,mesh.nfes,mesh.nfemax,mesh.ne);                 
//     }
//     else {
//         Int nmax = mesh.npfmax*mesh.nfemax;    
//         vector<Int> nmaxs = mesh.nfes;
//         Int nfe, nelem = mesh.nfes.size();
//         for(Int i=0;i<nelem;i++) {
//             nfe = mesh.nfes[i];
//             nmaxs[i] = 0;
//             for (Int j=0; j<nfe; j++)
//                 nmaxs[i] += mesh.npfs[i][j];                 
//         }
//         elcon2entcon(rowent2elem,colent2elem,sys.ent2entStart,sys.ent2ent,mesh.elcon,mesh.elementtype,nmaxs,nmax,mesh.ne);         
//     }
//     
//     sys.nproc = app.nproc;
//     if (hybrid == 0) //HDG
//         sys.blkSize = app.nch*mesh.npfmax;
//     else
//         sys.blkSize = app.nch;    
//     sys.BJ_nrows = sys.ent2entStart.size()-1;
//     sys.numEntities = sys.ent2entStart.size()-1;
//     sys.numBlocks = sys.ent2ent.size();    
//         
//     Int a;
//     sys.maxBlocksPerRow = 0;
//     for (Int i=1; i<sys.ent2entStart.size(); i++) {
//         a = sys.ent2entStart[i]-sys.ent2entStart[i-1];
//         if (a > sys.maxBlocksPerRow)
//             sys.maxBlocksPerRow = a;
//     }    
// }
// 
// 
// #ifdef HAVE_MPI
// 
// void dmdparallel(sysstruct &sys, meshstruct &mesh, appstruct &app, solstruct &sol, dmdstruct &dmd)
// {
//     Int N = mesh.ndims.size();
//     Int nie = mesh.ndims[17];  // number of interior elements
//     Int nif = mesh.ndims[18];  // number of interior faces                    
//     Int nin = mesh.ndims[19];  // number of interior entities                     
//     Int nxe = mesh.ndims[21];  // number of exterior elements
//     Int nxf = mesh.ndims[22];  // number of exterior faces     
//     Int nxn = mesh.ndims[23];  // number of exterior entities            
//     Int hybrid = mesh.ndims[N-2];                
//     
//     if (hybrid == 0) {  // HDG              
//         domaindecomposition(dmd, mesh.t2f, mesh.elcon, mesh.t2f, mesh.elementtype, mesh.extintelem, mesh.extintface,
//                 mesh.elem2cpu, mesh.face2cpu, mesh.nfes, mesh.npfs, nie, nin, nxe, nxn, 1, mesh.my_rank);                        
//         Int nnbsd = (Int) dmd.nbsd.size();                  
//         vector< dmdstruct > nbdmd(nnbsd, dmdstruct());
//         make_nbdmd(nbdmd, dmd);        
//         sendrecvhdg(dmd);
//     }
//     else {                        
//         domaindecomposition(dmd, mesh.elcon, mesh.elcon, mesh.t2f, mesh.elementtype, mesh.extintelem, mesh.extintent,
//                 mesh.elem2cpu, mesh.ent2cpu, mesh.nfes, mesh.npfs, nie, nin, nxe, nxn, 0, mesh.my_rank);                
//         Int nnbsd = (Int) dmd.nbsd.size();          
//         vector< dmdstruct > nbdmd(nnbsd, dmdstruct());
//         make_nbdmd(nbdmd, dmd);        
//         sendrecvedg(dmd, nbdmd);
//     }                        
//     
//     // update mesh
//     mesh.ne = dmd.elemmap.size();    
//     mesh.nf = *max_element(dmd.t2f.begin(), dmd.t2f.end())+1;
//     mesh.ndofuh = *max_element(dmd.elcon.begin(), dmd.elcon.end())+1;
//     
//     mesh.t2f = dmd.t2f; //iarray2datn(dmd.t2f, mesh.ne, mesh.nfemax);
//     mesh.elcon = dmd.elcon; //iarray2datn(dmd.elcon, ne, mesh.nfemax*mesh.npfmax);
// 
//     vector<double> tm; 
//     vector<Int> tn;
//     tm = darray2datindex(sol.UDG, dmd.elemmap, mesh.npemax*app.nc); sol.UDG = tm;
//     tm = darray2datindex(mesh.dgnodes, dmd.elemmap, mesh.nmemax*app.ncd); mesh.dgnodes = tm;           
//     tn = iarray2datindex(mesh.bf, dmd.elemmap, mesh.nfemax); mesh.bf = tn;    
//     tn = iarrayatindex(mesh.elementtype, dmd.elemmap); mesh.elementtype = tn;
//     tn = iarray2datindex(mesh.t, dmd.elemmap, mesh.nvemax); mesh.t = tn;  
//     tn = iarray2datindex(mesh.t2t, dmd.elemmap, mesh.nfemax); mesh.t2t = tn;   
//     tn = iarrayatindex(mesh.elem2cpu, dmd.elemmap); mesh.elem2cpu = tn;
//     tn = iarrayatindex(mesh.extintelem, dmd.elemmap); mesh.extintelem = tn;   
//     
// //     if (sys.my_rank==0) 
// //         print1iarray(&mesh.bf[0],2*mesh.nfemax);
//         
//     if (hybrid == 0)  { // HDG 
// //         if (app.my_rank==0) {
// //             print2darray(&sol.UH[0], app.nch, mesh.npfmax);
// //         }
//         tm = darray2datindex(sol.UH, dmd.entmap, app.nch*mesh.npfmax); 
//         sol.UH = tm;   
// //         if (app.my_rank==0) {
// //             print2darray(&sol.UH[0], app.nch, mesh.npfmax);
// //         }
// //         error("here");
//     }
//     else {
//         tm = darray2datindex(sol.UH, dmd.entmap, app.nch); 
//         sol.UH = tm;         
//     }
//     
//     // update sys    
//     sys.entpart = dmd.entpart;
//     sys.entpartpts = dmd.entpartpts;
//     sys.entrecv = dmd.entrecv;
//     sys.entrecvpts = dmd.entrecvpts;
//     sys.entsend = dmd.entsend;
//     sys.entsendpts = dmd.entsendpts;
//     sys.elempart = dmd.elempart;
//     sys.elempartpts = dmd.elempartpts;
//     sys.elemrecv = dmd.elemrecv;
//     sys.elemrecvpts = dmd.elemrecvpts;
//     sys.elemsend = dmd.elemsend;
//     sys.elemsendpts = dmd.elemsendpts;
//     sys.vecrecv = dmd.vecrecv;
//     sys.vecrecvpts = dmd.vecrecvpts;
//     sys.vecsend = dmd.vecsend;
//     sys.vecsendpts = dmd.vecsendpts;    
//     sys.matrecv = dmd.matrecv;
//     sys.matrecvpts = dmd.matrecvpts;
//     sys.matsend = dmd.matsend;
//     sys.matsendpts = dmd.matsendpts;    
//     sys.ent2entStart = dmd.bcrs_rowent2ent;
//     sys.ent2ent = dmd.bcrs_colent2ent;
//     sys.nbsd = dmd.nbsd;
//     
//     sys.nproc = app.nproc;
//     sys.my_rank = app.my_rank;    
//     if (hybrid == 0) //HDG
//         sys.blkSize = app.nch*mesh.npfmax;
//     else
//         sys.blkSize = app.nch;    
//     sys.nentrecv = sys.entrecv.size();
//     sys.nentsend = sys.entsend.size();    
//     sys.nelemrecv = sys.elemrecv.size();
//     sys.nelemsend = sys.elemsend.size();
//     sys.nvecrecv = sys.matrecv.size();
//     sys.nvecsend = sys.matsend.size();
//     sys.nmatrecv = sys.matrecv.size();
//     sys.nmatsend = sys.matsend.size();    
//     sys.nnbsd = sys.nbsd.size();    
//     sys.BJ_nrows = sys.entpartpts[0]+sys.entpartpts[1];
//     sys.numEntities = sys.entpart.size();
//     sys.numBlocks = sys.ent2ent.size();    
//     
//     Int a;
//     sys.maxBlocksPerRow = 0;
//     for (Int i=1; i<sys.ent2entStart.size(); i++) {
//         a = sys.ent2entStart[i]-sys.ent2entStart[i-1];
//         if (a > sys.maxBlocksPerRow)
//             sys.maxBlocksPerRow = a;
//     }
// }    
// 
// #endif


//
// vector<Int> getndims(meshstruct &mesh, vector< masterstruct > &master, appstruct &app, solstruct &sol, sysstruct &sys)
// {
//     vector<Int> ndims(42,0);
//     ndims[0] = app.nd;
//     ndims[1] = app.ncd;
//     ndims[2] = app.nc;
//     ndims[3] = app.ncu;
//     ndims[4] = app.ncq;
//     ndims[5] = app.ncp;
//     ndims[6] = app.nch;
//     ndims[7] = app.nco;    
//     ndims[8] = mesh.ne;
//     ndims[9] = mesh.nf;
//     ndims[10] = mesh.nv;
//     ndims[11] = mesh.ndofuh;
//     ndims[12] = mesh.nfemax;
//     ndims[13] = mesh.nvemax;
//     ndims[14] = mesh.nlemax;
//     ndims[15] = mesh.npemax;
//     ndims[16] = mesh.nmemax;
//     ndims[17] = mesh.ngemax;
//     ndims[18] = mesh.ngeRmax;
//     ndims[19] = mesh.nvfmax;    
//     ndims[20] = mesh.npfmax;    
//     ndims[21] = mesh.nmfmax;    
//     ndims[22] = mesh.ngfmax;
//     ndims[23] = mesh.ngfRmax;
//     ndims[24] = mesh.ndfmax;
//     ndims[25] = mesh.ncfmax;
//     ndims[26] = sys.numEntities;
//     ndims[27] = sys.numBlocks;
//     ndims[28] = sys.maxBlocksPerRow;
//     ndims[29] = sys.blkSize;
//     ndims[30] = sys.BJ_nrows;    
//     ndims[31] = sys.nproc;    
//     ndims[32] = sys.nnbsd;
//     ndims[33] = sys.nentrecv;
//     ndims[34] = sys.nentsend;    
//     ndims[35] = sys.nelemrecv;
//     ndims[36] = sys.nelemsend;    
//     ndims[37] = sys.nvecrecv;
//     ndims[38] = sys.nvecsend;                
//     ndims[39] = sys.nmatrecv;
//     ndims[40] = sys.nmatsend;                
//     ndims[41] = sys.my_rank;    
//     
//     return ndims;
// }

#endif