#ifndef __SENDRECVMPI
#define __SENDRECVMPI

// Written by: C. Nguyen & P. Fernandez

#ifdef  HAVE_MPI

void sendvector(sysstruct &sys, double* x, Int * requestCounter)
{
    Int inc = 1;
    
    // Initialize request counter:
    Int requestCounter_tmp = 0;
    
    // Copy some portion of x to buffsend:
    Int buff_pt = 0;
    for (Int i=0; i<sys.nentsend; i++) {
        Int ri = sys.entsend[i];
        Int x_pt = sys.rowStartR[ri];
        DCOPY(&entSz[ri], &x[x_pt], &inc, &sys.buffsend[buff_pt], &inc);
        buff_pt += entSz[ri];
    }
    
    // Non-blocking send:
    Int psend = 0;
    for (Int n=0; n<sys.nnbsd; n++) {
        Int neighbor = sys.nbsd[n];
        Int nsend = sys.entsendlen[n];
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            psend += nsend;
            requestCounter_tmp += 1;
        }
    }

    // Non-blocking receive:
    Int precv = 0;
    for (Int n=0; n<sys.nnbsd; n++) {
        Int neighbor = sys.nbsd[n];
        Int nrecv = sys.entrecvlen[n];
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            precv += nrecv;
            requestCounter_tmp += 1;
        }
    }

    *requestCounter = requestCounter_tmp;
}

void recvvector(sysstruct &sys, double* x, Int * requestCounter)
{
    Int inc = 1;
    
    // Wait until all send and receive operations are completely done:
    MPI_Waitall(*requestCounter, sys.requests, sys.statuses);
    
    // Copy buffrecv to x:
    Int buff_pt = 0;
    for (Int i=0; i<sys.nentrecv; i++) {
        Int ri = sys.entrecv[i];
        Int x_pt = rowStartR[ri];
        DCOPY(&entSz[ri], &sys.buffrecv[buff_pt], &inc, &x[x_pt], &inc);
        buff_pt += entSz[ri];
    }
}

void sendrecvvector(sysstruct &sys, double* x)
{            
    Int requestCounter[1];
    sendvector(sys, &x[0], &requestCounter[0]);
    recvvector(sys, &x[0], &requestCounter[0]);
}

void sendUDG(sysstruct &sys, double* UDG, Int* ndims, Int* requestCounter)
{
    // UDG: npe / nc / ne
    
    // Initialize request counter:
    Int requestCounter_tmp = 0;
    
    // Copy some portion of UDG to buffsend:
    Int buff_pt = 0;
    for (Int i=0; i<sys.nelemsend; i++) {
        Int ie = sys.elemsend[i];
        Int elementype = mesh.elementype[ie];
        Int npv = mesh.elements[elementype].npv;
        Int nc = ndims[19];
        Int len = npv*nc;
        Int UDG_pt = sys.UDGstart[ie];
        DCOPY(&len, &UDG[UDG_pt], &inc, &sys.buffsend[buff_pt], &inc);
        buff_pt += len;
    }

    // Non-blocking send:
    Int psend = 0;
    for (Int n=0; n<sys.nnbsd; n++) {
        Int neighbor = sys.nbsd[n];
        Int nsend = sys.elemsendlenUDG[n];
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            psend += nsend;
            requestCounter_tmp += 1;
        }
    }
    
    // Non-blocking receive:
    Int precv = 0;
    for (Int n=0; n<sys.nnbsd; n++) {
        Int neighbor = sys.nbsd[n];
        Int nrecv = sys.elemrecvlenUDG[n];
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            precv += nrecv;
            requestCounter_tmp += 1;
        }
    }

    *requestCounter = requestCounter_tmp;
}

void recvUDG(sysstruct &sys, double* UDG, Int* ndims, Int * requestCounter)
{
    // UDG: npe / nc / ne
    
    // Wait until all send and receive operations are completely done:
    MPI_Waitall(*requestCounter, sys.requests, sys.statuses);
    
    // Copy buffrecv to UDG:
    Int buff_pt = 0;
    for (Int i=0; i<sys.nelemrecv; i++) {
        Int ie = sys.elemrecv[i];
        Int elementype = mesh.elementype[ie];
        Int npv = mesh.elements[elementype].npv;
        Int nc = ndims[19];
        Int len = npv*nc;
        Int UDG_pt = sys.UDGstart[ie];
        DCOPY(&len, &sys.buffrecv[buff_pt], &inc, &UDG[UDG_pt], &inc);
        buff_pt += len;
    }
}

void sendrecvUDG(sysstruct &sys, double* UDG, Int* ndims)
{            
    Int requestCounter[1];
    sendUDG(sys, &UDG[0], &ndims[0], &requestCounter[0]);
    recvUDG(sys, &UDG[0], &ndims[0], &requestCounter[0]);
}

void sendU(sysstruct &sys, double* U, Int* ndims, Int* requestCounter)
{
    // U: npe / ncu / ne
    
    // Initialize request counter:
    Int requestCounter_tmp = 0;
    
    // Copy some portion of U to buffsend:
    Int buff_pt = 0;
    for (Int i=0; i<sys.nelemsend; i++) {
        Int ie = sys.elemsend[i];
        Int elementype = mesh.elementype[ie];
        Int npv = mesh.elements[elementype].npv;
        Int ncu = ndims[20];
        Int len = npv*ncu;
        Int U_pt = sys.Ustart[ie];
        DCOPY(&len, &U[U_pt], &inc, &sys.buffsend[buff_pt], &inc);
        buff_pt += len;
    }

    // Non-blocking send:
    Int psend = 0;
    for (Int n=0; n<sys.nnbsd; n++) {
        Int neighbor = sys.nbsd[n];
        Int nsend = sys.elemsendU[n];
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            psend += nsend;
            requestCounter_tmp += 1;
        }
    }
    
    // Non-blocking receive:
    Int precv = 0;
    for (Int n=0; n<sys.nnbsd; n++) {
        Int neighbor = sys.nbsd[n];
        Int nrecv = sys.elemrecvlenU[n];
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            precv += nrecv;
            requestCounter_tmp += 1;
        }
    }

    *requestCounter = requestCounter_tmp;
}

void recvU(sysstruct &sys, double* U, Int* ndims, Int * requestCounter)
{
    // U: npe / ncu / ne
    
    // Wait until all send and receive operations are completely done:
    MPI_Waitall(*requestCounter, sys.requests, sys.statuses);
    
    // Copy buffrecv to U:
    Int buff_pt = 0;
    for (Int i=0; i<sys.nelemrecv; i++) {
        Int ie = sys.elemrecv[i];
        Int elementype = mesh.elementype[ie];
        Int npv = mesh.elements[elementype].npv;
        Int ncu = ndims[20];
        Int len = npv*ncu;
        Int U_pt = sys.Ustart[ie];
        DCOPY(&len, &sys.buffrecv[buff_pt], &inc, &U[U_pt], &inc);
        buff_pt += len;
    }
}


void sendrecvU(sysstruct &sys, double* U, Int* ndims)
{        
    Int requestCounter[1];
    sendU(sys, &U[0], &ndims[0], &requestCounter[0]);
    recvU(sys, &U[0], &ndims[0], &requestCounter[0]);
}

void sendScalarField(sysstruct &sys, double* field, Int* ndims, Int* requestCounter)
{
    // field: npe / ne
    
    // Initialize request counter:
    Int requestCounter_tmp = 0;
    
    // Copy some portion of field to buffsend:
    Int buff_pt = 0;
    for (Int i=0; i<sys.nelemsend; i++) {
        Int ie = sys.elemsend[i];
        Int elementype = mesh.elementype[ie];
        Int npv = mesh.elements[elementype].npv;
        Int len = npv;
        Int field_pt = sys.scalarFieldStart[ie];
        DCOPY(&len, &field[field_pt], &inc, &sys.buffsend[buff_pt], &inc);
        buff_pt += len;
    }

    // Non-blocking send:
    Int psend = 0;
    for (Int n=0; n<sys.nnbsd; n++) {
        Int neighbor = sys.nbsd[n];
        Int nsend = sys.elemsendScalar[n];
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            psend += nsend;
            requestCounter_tmp += 1;
        }
    }
    
    // Non-blocking receive:
    Int precv = 0;
    for (Int n=0; n<sys.nnbsd; n++) {
        Int neighbor = sys.nbsd[n];
        Int nrecv = sys.elemrecvlenScalar[n];
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            precv += nrecv;
            requestCounter_tmp += 1;
        }
    }

    *requestCounter = requestCounter_tmp;
}

void recvScalarField(sysstruct &sys, double* field, Int* ndims, Int * requestCounter)
{
    // field: npe / ne
    
    // Wait until all send and receive operations are completely done:
    MPI_Waitall(*requestCounter, sys.requests, sys.statuses);
    
    // Copy buffrecv to field:
    Int buff_pt = 0;
    for (Int i=0; i<sys.nelemrecv; i++) {
        Int ie = sys.elemrecv[i];
        Int elementype = mesh.elementype[ie];
        Int npv = mesh.elements[elementype].npv;
        Int len = npv;
        Int field_pt = sys.scalarFieldStart[ie];
        DCOPY(&len, &sys.buffrecv[buff_pt], &inc, &field[field_pt], &inc);
        buff_pt += len;
    }
}

void sendrecvScalarField(sysstruct &sys, double* field, Int* ndims)
{
    Int requestCounter[1];
    sendScalarField(sys, &field[0], &ndims[0], &requestCounter[0]);
    recvScalarField(sys, &field[0], &ndims[0], &requestCounter[0]);
}

// // // // void sendrecvmatrix(sysstruct &sys, double* Hg)
// // // // {
// // // //     Int inc = 1;
// // // //     
// // // //     // Initialize request counter:
// // // //     Int request_counter = 0;
// // // //     
// // // //     // Copy some portion of Hg to buffsendmat:
// // // //     for (Int i=0; i<sys.nmatsend; i++) {
// // // //         Int pi = sys.matsend[i];
// // // //         Int bsz2 = entSz[matsend_rowEnt[i]] * entSz[matsend_colEnt[i]];
// // // //         DCOPY(&bsz2, &Hg[blkStartJ[pi]], &inc, &buffsendmat[buffsendmatStart[i]], &inc);
// // // //     }
// // // //     
// // // //     // Non-blocking send:
// // // //     Int psend = 0;
// // // //     for (Int n=0; n<sys.nnbsd; n++) {
// // // //         Int neighbor = sys.nbsd[n];
// // // //         Int nsend = sys.matsendlen[n];
// // // //         if (nsend>0) {
// // // //             MPI_Isend(&buffsendmat[psend], nsend, MPI_DOUBLE, neighbor, 0,
// // // //                    MPI_COMM_WORLD, &requests[request_counter]);
// // // //             psend += nsend;
// // // //             request_counter += 1;
// // // //         }
// // // //     }
// // // // 
// // // //     // Non-blocking receive:
// // // //     Int precv = 0;
// // // //     for (Int n=0; n<sys.nnbsd; n++) {
// // // //         Int neighbor = sys.nbsd[n];
// // // //         Int nrecv = sys.matrecvlen[n];
// // // //         if (nrecv>0) {
// // // //             MPI_Irecv(&buffrecvmat[precv], nrecv, MPI_DOUBLE, neighbor, 0,
// // // //                    MPI_COMM_WORLD, &requests[request_counter]);
// // // //             precv += nrecv;
// // // //             request_counter += 1;
// // // //         }
// // // //     }
// // // //     
// // // //     // Wait until all send and receive operations are completely done:
// // // //     MPI_Waitall(request_counter, requests, statuses);
// // // // 
// // // //     // Copy buffrecvmat to Hg:
// // // //     for (Int i=0; i<sys.nmatrecv; i++) {
// // // //         Int pi = sys.matrecv[i];
// // // //         Int bsz2 = entSz[matrecv_rowEnt[i]] * entSz[matrecv_colEnt[i]];
// // // //         DCOPY(&bsz2, &buffrecvmat[buffrecvmatStart[i]], &inc, &sys.Hg[blockStartJ[pi]], &inc);
// // // //     }
// // // // }

#endif

#endif
