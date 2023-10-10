#ifndef __SENDRECVMPI
#define __SENDRECVMPI

// Written by: C. Nguyen & P. Fernandez

#ifdef  HAVE_MPI

static void sendvector(sysstruct &sys, double* x, Int * requestCounter)
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
            MPI_Isend(&sys.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
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
            MPI_Irecv(&sys.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                      MPI_COMM_WORLD, &sys.requests[requestCounter_tmp]);
            precv += nrecv;
            requestCounter_tmp += 1;
        }
    }

    *requestCounter = requestCounter_tmp;
}

static void recvvector(sysstruct &sys, double* x, Int * requestCounter)
{
    Int i, ri, inc = 1;
    Int bsz = sys.blkSize;

    /* Wait until all send and receive operations are completely done */
    MPI_Waitall(*requestCounter, sys.requests, sys.statuses);

    /* Copy buffrecv to x */
    for (i=0; i<sys.nentrecv; i++) {
        ri = sys.entrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &x[bsz*ri], &inc);
    }
}


static void sendrecvvector(sysstruct &sys, double* x)
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
            MPI_Isend(&sys.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
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
            MPI_Irecv(&sys.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[requestCounter]);
            precv += nrecv;
            requestCounter += 1;
        }        
    }    
    
    /* Wait until all send and receive operations are completely done */
    MPI_Waitall(requestCounter, sys.requests, sys.statuses);
        
    /* Copy buffrecv to x */
    for (i=0; i<sys.nentrecv; i++) {
        ri = sys.entrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &x[bsz*ri], &inc);
    }    
}

#endif

#endif
