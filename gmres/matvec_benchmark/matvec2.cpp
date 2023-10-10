#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

void matvec2(float* v, float *Hg, float* x, sysstruct &sys, int preftchDist)
{
//     DESCRIPTION: v = Hg*x, where Hg is a sparse matrix and (x,v) are vectors. Single precision implementation.
//     Implementation that performs the all the small matrix-vector products in a given row at once (one SGEMV call only)
    
//     INPUTS / OUTPUTS:
//     Hg:      Pointer to the sparse matrix Hg.
//              Hg is stored in Block Compressed Row Storage (BCRS) format. More info: http://www.netlib.org/utk/people/JackDongarra/etemplates/node375.html
//     x:       Pointer to the input vector x.
//     v:       Pointer to the output vector v.
//     sys:     "sysstruct" structure with information about the linear system to be solved.

//     NOTES:
//     - SGEMV is the function of the BLAS library to perform matrix-vector products in single precision. DGEMV is the analogous function for double precision.
//     - For our purpose here, we can assume sys.blkSize = 5

	float one = 1.0, zero = 0.0, alpha;
	char chn = 'N';

	Int nri, pij, pii, ri, rj, nb, inc = 1;
	Int bsz = sys.blkSize;              // Size of each block in Hg. Let's assume bsz = 5
	Int bsz2 = bsz*bsz;

	// Row block contains pointers to the start of each row.
	Int *row_blk = &sys.row_blk[0];     // Description of this vector is in: http://www.netlib.org/utk/people/JackDongarra/etemplates/node375.html
    
	// Col block contains column indices of each non-zero block
	Int *col_ind = &sys.col_ind[0];     // Description of this vector is in: http://www.netlib.org/utk/people/JackDongarra/etemplates/node375.html

	clock_t t0 = clock();

    for (Int ri = 0; ri < sys.numRows; ri++) {
        nri = row_blk[ri+1] - row_blk[ri];  // Number of non-zero blocks in current row
        pii = row_blk[ri];
        pij = row_blk[ri] - 1;
        
//         __builtin_prefetch(&Hg[bsz2*(pii+preftchDist)] , 0,3);
        
        // Create temporary dense vector (to perform matrix-vector product of the entire row with one BLAS call only):
        for (Int i = 0; i<nri; i++) {
            pij += 1;
            rj = col_ind[pij];
            for (Int j = 0; j < bsz; j++) {
//                 __builtin_prefetch(&x[bsz*col_ind[pij+preftchDist]], 0,3);
                sys.xDenseRow[bsz*i+j] = x[bsz*rj+j];
            }
        }
        
        // Perform matrix-vector product using BLAS function SGEMV:
        nb = bsz*nri;
        SGEMV(
            &chn, // type of transformation: y := alpha*A*x + beta*y.
            &bsz, // num rows
            &nb, // num columns
            &one, // scalar alpha, is 1.0
            &Hg[bsz2*pii], // start of current row
            &bsz, // dimension of block
            &sys.xDenseRow[0], // 
            &inc, // inc of lee n y
            &zero, // scalar beta, 0.0
            &v[bsz*ri], // location of y to write at
            &inc // inc of elements in y
        );
    }
	double tTotal = clock() - t0;
}
