#ifndef __SENDRECVMPI
#define __SENDRECVMPI

// Written by: C. Nguyen & P. Fernandez

#ifdef  HAVE_MPI

void sendvector(sysstruct &sys, double* x, Int * requestCounter)
{
    Int i, n, ri, neighbor, inc = 1;
    Int bsz = sys.blkSize;

    /* Initialize request counter */
    Int requestCounter_tmp = 0;

    /* Copy some portion of x to buffsend */
    for (i=0; i<sys.nentsend; i++) {
        ri = sys.entsend[i];
        DCOPY(&bsz, &x[bsz*ri], &inc, &sys.buffsend[bsz*i], &inc);
    }

    /* Non-blocking send */
    Int nsend, psend = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nsend = sys.entsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], (int) nsend, MPI_DOUBLE, (int) neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            psend += nsend;
            requestCounter_tmp += 1;
        }
    }

    /* Non-blocking receive */
    Int nrecv, precv = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nrecv = sys.entrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], (int) nrecv, MPI_DOUBLE, (int) neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            precv += nrecv;
            requestCounter_tmp += 1;
        }
    }

    *requestCounter = requestCounter_tmp;
}

void recvvector(sysstruct &sys, double* x, Int * requestCounter)
{
    Int i, ri, inc = 1;
    Int bsz = sys.blkSize;

    /* Wait until all send and receive operations are completely done */
    MPI_Waitall((int) (*requestCounter), sys.requests, sys.statuses);

    /* Copy buffrecv to x */
    for (i=0; i<sys.nentrecv; i++) {
        ri = sys.entrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &x[bsz*ri], &inc);
    }
}


void sendrecvvector(sysstruct &sys, double* x)
{            
    Int i, n, ri, neighbor, inc = 1;
    Int bsz = sys.blkSize;

    /* Initialize request counter */
    Int requestCounter = 0;        
    
    /* Copy some portion of x to buffsend */
    for (i=0; i<sys.nentsend; i++) {
        ri = sys.entsend[i];
        DCOPY(&bsz, &x[bsz*ri], &inc, &sys.buffsend[bsz*i], &inc);
    }

    /* Non-blocking send */
    Int nsend, psend = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nsend = sys.entsendpts[n]*bsz;        
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], (int) nsend, MPI_DOUBLE, (int) neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[requestCounter]);
            psend += nsend;
            requestCounter += 1;
        }        
    }    
    
    /* Non-blocking receive */
    Int nrecv, precv = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nrecv = sys.entrecvpts[n]*bsz;        
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], (int) nrecv, MPI_DOUBLE, (int) neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[requestCounter]);
            precv += nrecv;
            requestCounter += 1;
        }        
    }    
    
    /* Wait until all send and receive operations are completely done */
    MPI_Waitall((int) requestCounter, sys.requests, sys.statuses);
        
    /* Copy buffrecv to x */
    for (i=0; i<sys.nentrecv; i++) {
        ri = sys.entrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &x[bsz*ri], &inc);
    }    
}


void sendUDG(sysstruct &sys, double* UDG, Int* ndims, Int* requestCounter)
{
    Int i, n, ri, neighbor, inc = 1;
    Int npe = ndims[9];
    Int nc = ndims[19];
    Int bsz = npe*nc;

    /* Initialize request counter */
    Int requestCounter_tmp = 0;

    /* Copy some portion of UDG to buffsend */
    for (i=0; i<sys.nelemsend; i++) {
        ri = sys.elemsend[i];
        DCOPY(&bsz, &UDG[bsz*ri], &inc, &sys.buffsend[bsz*i], &inc);
    }

    /* Non-blocking send */
    Int nsend, psend = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nsend = sys.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], (int) nsend, MPI_DOUBLE, (int) neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            psend += nsend;
            requestCounter_tmp += 1;
        }
    }

    /* Non-blocking receive */
    Int nrecv, precv = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nrecv = sys.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], (int) nrecv, MPI_DOUBLE, (int) neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            precv += nrecv;
            requestCounter_tmp += 1;
        }
    }

    *requestCounter = requestCounter_tmp;
}


void recvUDG(sysstruct &sys, double* UDG, Int* ndims, Int * requestCounter)
{
    Int i, ri, inc = 1;
    Int npe = ndims[9];
    Int nc = ndims[19];
    Int bsz = npe*nc;

    /* Wait until all send and receive operations are completely done */
    MPI_Waitall((int) (*requestCounter), sys.requests, sys.statuses);

    /* Copy buffrecv to UDG */
    for (i=0; i<sys.nelemrecv; i++) {
        ri = sys.elemrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &UDG[bsz*ri], &inc);
    }
}


void sendrecvUDG(sysstruct &sys, double* UDG, Int* ndims)
{            
    Int i, n, ri, neighbor, inc = 1;
    Int npe = ndims[9];
    Int nc = ndims[19];
    Int bsz = npe*nc;

    /* Initialize request counter */
    Int requestCounter = 0;        
        
    /* Copy some portion of UDG to buffsend */
    for (i=0; i<sys.nelemsend; i++) {
        ri = sys.elemsend[i];
        DCOPY(&bsz, &UDG[bsz*ri], &inc, &sys.buffsend[bsz*i], &inc);
    }

    /* Non-blocking send */
    Int nsend, psend = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nsend = sys.elemsendpts[n]*bsz;        
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], (int) nsend, MPI_DOUBLE, (int) neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[requestCounter]);
            psend += nsend;
            requestCounter += 1;
        }        
    }    
    
    /* Non-blocking receive */
    Int nrecv, precv = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nrecv = sys.elemrecvpts[n]*bsz;        
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], (int) nrecv, MPI_DOUBLE, (int) neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[requestCounter]);
            precv += nrecv;
            requestCounter += 1;
        }        
    }    
    
    /* Wait until all send and receive operations are completely done */
    MPI_Waitall((int) requestCounter, sys.requests, sys.statuses);
    
    /* Copy buffrecv to UDG */
    for (i=0; i<sys.nelemrecv; i++) {
        ri = sys.elemrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &UDG[bsz*ri], &inc);
    }    
}


void sendU(sysstruct &sys, double* U, Int* ndims, Int* requestCounter)
{
    // U: npe / ncu / ne
    
    Int i, n, ri, neighbor, inc = 1;
    Int npe = ndims[9];
    Int ncu = ndims[20];
    Int bsz = npe*ncu;
    
    /* Initialize request counter */
    Int requestCounter_tmp = 0;

    /* Copy some portion of U to buffsend */
    for (i=0; i<sys.nelemsend; i++) {
        ri = sys.elemsend[i];
        DCOPY(&bsz, &U[bsz*ri], &inc, &sys.buffsend[bsz*i], &inc);
    }

    /* Non-blocking send */
    Int nsend, psend = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nsend = sys.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], (int) nsend, MPI_DOUBLE, (int) neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            psend += nsend;
            requestCounter_tmp += 1;
        }
    }

    /* Non-blocking receive */
    Int nrecv, precv = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nrecv = sys.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], (int) nrecv, MPI_DOUBLE, (int) neighbor, 0,
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
    
    Int i, ri, inc = 1;
    Int npe = ndims[9];
    Int ncu = ndims[20];
    Int bsz = npe*ncu;

    /* Wait until all send and receive operations are completely done */
    MPI_Waitall((int) (*requestCounter), sys.requests, sys.statuses);

    /* Copy buffrecv to U */
    for (i=0; i<sys.nelemrecv; i++) {
        ri = sys.elemrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &U[bsz*ri], &inc);
    }
}


void sendrecvU(sysstruct &sys, double* U, Int* ndims)
{           
    // U: npe / ncu / ne
    
    Int i, n, ri, neighbor, inc = 1;
    Int npe = ndims[9];
    Int ncu = ndims[20];
    Int bsz = npe*ncu;

    /* Initialize request counter */
    Int requestCounter = 0;        
        
    /* Copy some portion of U to buffsend */
    for (i=0; i<sys.nelemsend; i++) {
        ri = sys.elemsend[i];
        DCOPY(&bsz, &U[bsz*ri], &inc, &sys.buffsend[bsz*i], &inc);
    }

    /* Non-blocking send */
    Int nsend, psend = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nsend = sys.elemsendpts[n]*bsz;        
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], (int) nsend, MPI_DOUBLE, (int) neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[requestCounter]);
            psend += nsend;
            requestCounter += 1;
        }        
    }    
    
    /* Non-blocking receive */
    Int nrecv, precv = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nrecv = sys.elemrecvpts[n]*bsz;        
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], (int) nrecv, MPI_DOUBLE, (int) neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[requestCounter]);
            precv += nrecv;
            requestCounter += 1;
        }        
    }    
    
    /* Wait until all send and receive operations are completely done */
    MPI_Waitall((int) requestCounter, sys.requests, sys.statuses);
    
    /* Copy buffrecv to U */
    for (i=0; i<sys.nelemrecv; i++) {
        ri = sys.elemrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &U[bsz*ri], &inc);
    }    
}

void sendScalarField(sysstruct &sys, double* field, Int* ndims, Int* requestCounter)
{
    // field: npe / ne
    
    Int i, n, ri, neighbor, inc = 1;
    Int npe = ndims[9];
    Int bsz = npe;
    
    /* Initialize request counter */
    Int requestCounter_tmp = 0;

    /* Copy some portion of field to buffsend */
    for (i=0; i<sys.nelemsend; i++) {
        ri = sys.elemsend[i];
        DCOPY(&bsz, &field[bsz*ri], &inc, &sys.buffsend[bsz*i], &inc);
    }

    /* Non-blocking send */
    Int nsend, psend = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nsend = sys.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], (int) nsend, MPI_DOUBLE, (int) neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            psend += nsend;
            requestCounter_tmp += 1;
        }
    }

    /* Non-blocking receive */
    Int nrecv, precv = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nrecv = sys.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], (int) nrecv, MPI_DOUBLE, (int) neighbor, 0,
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
    
    Int i, ri, inc = 1;
    Int npe = ndims[9];
    Int bsz = npe;

    /* Wait until all send and receive operations are completely done */
    MPI_Waitall((int) (*requestCounter), sys.requests, sys.statuses);

    /* Copy buffrecv to field */
    for (i=0; i<sys.nelemrecv; i++) {
        ri = sys.elemrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &field[bsz*ri], &inc);
    }
}

void sendrecvScalarField(sysstruct &sys, double* field, Int* ndims)
{           
    // field: npe / ne
    
    Int i, n, ri, neighbor, inc = 1;
    Int npe = ndims[9];
    Int bsz = npe;

    /* Initialize request counter */
    Int requestCounter = 0;        
        
    /* Copy some portion of field to buffsend */
    for (i=0; i<sys.nelemsend; i++) {
        ri = sys.elemsend[i];
        DCOPY(&bsz, &field[bsz*ri], &inc, &sys.buffsend[bsz*i], &inc);
    }

    /* Non-blocking send */
    Int nsend, psend = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nsend = sys.elemsendpts[n]*bsz;        
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], (int) nsend, MPI_DOUBLE, (int) neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[requestCounter]);
            psend += nsend;
            requestCounter += 1;
        }        
    }    
    
    /* Non-blocking receive */
    Int nrecv, precv = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nrecv = sys.elemrecvpts[n]*bsz;        
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], (int) nrecv, MPI_DOUBLE, (int) neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[requestCounter]);
            precv += nrecv;
            requestCounter += 1;
        }        
    }    
    
    /* Wait until all send and receive operations are completely done */
    MPI_Waitall((int) requestCounter, sys.requests, sys.statuses);
    
    /* Copy buffrecv to field */
    for (i=0; i<sys.nelemrecv; i++) {
        ri = sys.elemrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &field[bsz*ri], &inc);
    }    
}

void sendrecvmatrix(sysstruct &sys, double* Hg)
{
    Int i, n, ri, pi, neighbor, inc = 1;
    Int bsz = sys.blkSize*sys.blkSize;

    /* Initialize request counter */
    Int request_counter = 0;

    /* Copy some portion of Hg to buffsendmat */
    for (i=0; i<sys.nentsend; i++) {
        ri = sys.entsend[i];
        pi = sys.ent2entStart[ri];
        DCOPY(&bsz, &Hg[bsz*pi], &inc, &sys.buffsendmat[bsz*i], &inc);
    }

    /* Non-blocking send */
    Int nsend, psend = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nsend = sys.entsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&sys.buffsendmat[psend], (int) nsend, MPI_DOUBLE, (int) neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[request_counter]);
            psend += nsend;
            request_counter += 1;
        }        
    }
    
    /* Non-blocking receive */
    Int nrecv, precv = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nrecv = sys.entrecvpts[n]*bsz;        
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecvmat[precv], (int) nrecv, MPI_DOUBLE, (int) neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }        
    }    
    
    /* Wait until all send and receive operations are completely done */
    MPI_Waitall((int) request_counter, sys.requests, sys.statuses);
        
    /* Copy buffrecvmat to Hg */
    for (i=0; i<sys.nentrecv; i++) {
        ri = sys.entrecv[i];
        pi = sys.ent2entStart[ri];
        DCOPY(&bsz, &sys.buffrecvmat[bsz*i], &inc, &Hg[bsz*pi], &inc);
    }
}

#endif

#endif
