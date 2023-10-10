#ifndef __DMDPARLIB
#define __DMDPARLIB

// Written by: C. Nguyen & P. Fernandez

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <mpi.h>

#ifndef HAVE_MPI
#define HAVE_MPI
#endif

using namespace std;

#include "digasopre.h"
#include "errormsg.cpp"
#include "dmdpar.cpp"

// #include "mkdmd.cpp"
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
//     tn = iarrayatindex(mesh.isEDGelement, dmd.elemmap); mesh.isEDGelement = tn;
//     tn = iarray2datindex(mesh.t, dmd.elemmap, mesh.nvemax); mesh.t = tn;  
//     tn = iarray2datindex(mesh.t2t, dmd.elemmap, mesh.nfemax); mesh.t2t = tn;   
//     tn = iarrayatindex(mesh.elem2cpu, dmd.elemmap); mesh.elem2cpu = tn;
//     tn = iarrayatindex(mesh.extintelem, dmd.elemmap); mesh.extintelem = tn;   
//             
//     if (hybrid == 0)  { // HDG 
//         tm = darray2datindex(sol.UH, dmd.entmap, app.nch*mesh.npfmax); 
//         sol.UH = tm;   
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

// void createdmdstruct(dmdstruct &dmd, meshstruct &mesh)
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
// }    

// void setupsys(sysstruct &sys, meshstruct &mesh, appstruct &app, solstruct &sol)
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
//     dmdstruct dmd;
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
// //     mesh.dgnodes = darray2datindex(mesh.dgnodes, dmd.elemmap, mesh.nmemax*app.ncd);            
// //     mesh.bf = iarray2datindex(mesh.bf, dmd.elemmap, mesh.nfemax);     
// //     mesh.elementtype = iarrayatindex(mesh.elementtype, dmd.elemmap);
// //     mesh.t = iarray2datindex(mesh.t, dmd.elemmap, mesh.nvemax);   
// //     mesh.t2t = iarray2datindex(mesh.t2t, dmd.elemmap, mesh.nfemax);   
// //     mesh.elem2cpu = iarrayatindex(mesh.elem2cpu, dmd.elemmap);
// //     mesh.extintelem = iarrayatindex(mesh.extintelem, dmd.elemmap);    
// //     
// //     // update sol
// //     sol.UDG = darray2datindex(sol.UDG, dmd.elemmap, mesh.npemax*app.nc);         
// //     if (hybrid == 0)   // HDG 
// //         sol.UH = darray2datindex(sol.UH, dmd.entmap, mesh.npfmax*app.nch);         
// //     else
// //         sol.UH = darray2datindex(sol.UH, dmd.entmap, app.nch);         
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

// static void writeiarray(ofstream &out, vector<Int> &a)
// {
//     Int N = (Int) a.size();
//     if (N > 0)
//         out.write( reinterpret_cast<char*>( &a[0] ), sizeof(Int) * N );
// }
// 
// void writedmdstruct(string filename, dmdstruct &dmd)
// {
//     // Open file to write
//     ofstream out(filename.c_str(), ios::out | ios::binary);
// 
//     if (!out) {
//         error("Unable to open file " + filename);
//     }
// 
//     if (out) {                                
//         Int N = 50; 
//         vector<Int> ndims(N,0);                    
//         ndims[0] = N;
//         ndims[1] = dmd.my_rank;
//         ndims[2] = (Int) dmd.nbsd.size();
//         ndims[3] = (Int) dmd.intelem.size();
//         ndims[4] = (Int) dmd.intent.size();
//         ndims[5] = (Int) dmd.elempart.size();
//         ndims[6] = (Int) dmd.elempartpts.size();
//         ndims[7] = (Int) dmd.entpart.size();
//         ndims[8] = (Int) dmd.entpartpts.size();
//         ndims[9] = (Int) dmd.elemrecv.size();
//         ndims[10] = (Int) dmd.elemrecvpts.size();
//         ndims[11] = (Int) dmd.elemsend.size();
//         ndims[12] = (Int) dmd.elemsendpts.size();        
//         ndims[13] = (Int) dmd.entrecv.size();
//         ndims[14] = (Int) dmd.entrecvpts.size();
//         ndims[15] = (Int) dmd.entsend.size();
//         ndims[16] = (Int) dmd.entsendpts.size();
//         ndims[17] = (Int) dmd.vecrecv.size();
//         ndims[18] = (Int) dmd.vecrecvpts.size();
//         ndims[19] = (Int) dmd.vecsend.size();
//         ndims[20] = (Int) dmd.vecsendpts.size();
//         ndims[21] = (Int) dmd.matrecv.size();
//         ndims[22] = (Int) dmd.matrecvpts.size();
//         ndims[23] = (Int) dmd.matsend.size();
//         ndims[24] = (Int) dmd.matsendpts.size();
//         ndims[25] = (Int) dmd.rowent2elem.size();
//         ndims[26] = (Int) dmd.colent2elem.size();
//         ndims[27] = (Int) dmd.rowent2ent.size();
//         ndims[28] = (Int) dmd.colent2ent.size();
//         ndims[29] = (Int) dmd.bcrs_rowent2elem.size();
//         ndims[30] = (Int) dmd.bcrs_colent2elem.size();
//         ndims[31] = (Int) dmd.bcrs_rowent2ent.size();
//         ndims[32] = (Int) dmd.bcrs_colent2ent.size();        
//         ndims[33] = (Int) dmd.ent2ind.size();
//         ndims[34] = (Int) dmd.elcon.size();
//         ndims[35] = (Int) dmd.t2f.size();
//         ndims[36] = (Int) dmd.elemmap.size();
//         ndims[37] = (Int) dmd.entmap.size();
//         
//         /* Write dmd structure to a file */    
//         writeiarray(out, ndims);
//         writeiarray(out, dmd.nbsd);
//         writeiarray(out, dmd.intelem);
//         writeiarray(out, dmd.intent);
//         writeiarray(out, dmd.elempart);
//         writeiarray(out, dmd.elempartpts);
//         writeiarray(out, dmd.entpart);
//         writeiarray(out, dmd.entpartpts);
//         writeiarray(out, dmd.elemrecv);
//         writeiarray(out, dmd.elemrecvpts);
//         writeiarray(out, dmd.elemsend);
//         writeiarray(out, dmd.elemsendpts);
//         writeiarray(out, dmd.entrecv);
//         writeiarray(out, dmd.entrecvpts);
//         writeiarray(out, dmd.entsend);
//         writeiarray(out, dmd.entsendpts);
//         writeiarray(out, dmd.vecrecv);
//         writeiarray(out, dmd.vecrecvpts);
//         writeiarray(out, dmd.vecsend);
//         writeiarray(out, dmd.vecsendpts);
//         writeiarray(out, dmd.matrecv);
//         writeiarray(out, dmd.matrecvpts);
//         writeiarray(out, dmd.matsend);
//         writeiarray(out, dmd.matsendpts);
//         writeiarray(out, dmd.rowent2elem);
//         writeiarray(out, dmd.colent2elem);
//         writeiarray(out, dmd.rowent2ent);
//         writeiarray(out, dmd.colent2ent);
//         writeiarray(out, dmd.bcrs_rowent2elem);
//         writeiarray(out, dmd.bcrs_colent2elem);
//         writeiarray(out, dmd.bcrs_rowent2ent);
//         writeiarray(out, dmd.bcrs_colent2ent);
//         writeiarray(out, dmd.ent2ind);
//         writeiarray(out, dmd.elcon);
//         writeiarray(out, dmd.t2f);
//         writeiarray(out, dmd.elemmap);
//         writeiarray(out, dmd.entmap);
//     }
// 
//     // Close file:
//     out.close();
// }

#endif
