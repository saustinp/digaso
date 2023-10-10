#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

void matvec1(float* v, float *Hg, float* x, sysstruct &sys)
{
//     DESCRIPTION: v = Hg*x, where Hg is a sparse matrix and (x,v) are vectors. Single precision implementation.
//      Implementation that performs each small matrix-vector product individually.

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
    for (Int ri = 0; ri < sys.numRows; ri++) { // iterate over the rows
        alpha = 0.0;
        nri = row_blk[ri+1] - row_blk[ri];  // Number of non-zero blocks in current row
        // printf("%d: %d - %d\n", ri, row_blk[ri + 1], row_blk[ri]);
        for (Int i = 0; i<nri; i++) {
            pij = row_blk[ri]+i; // number of current block
            rj = col_ind[pij]; // column index of (1, 1) current block

            // Perform matrix-vector product using BLAS function SGEMV:
            // printf("hg: %d, x: %d, v: %d\n", bsz2 * pij, bsz * rj, bsz * ri);
            SGEMV(
                &chn, // type of transformation: y := alpha*A*x + beta*y.
                &bsz, // num rows
                &bsz, // num columns
                &one, // scalar alpha, is 1.0
                &Hg[bsz2*pij], // start of current block
                &bsz, // dimension of block
                &x[bsz*rj], // location of x to multiply into
                &inc, // increment for elements of x, is 1
                &alpha, // scalar beta, is 0 but changes to 1 to keep updating
                &v[bsz*ri], // location in y for output
                &inc // inc of elements in y, is 1
            );
            alpha = 1.0;
        }
    }
	double tTotal = clock() - t0;
}
