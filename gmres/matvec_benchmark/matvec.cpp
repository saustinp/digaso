#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

#include "wrapper.h"

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

void matvec2(float* v, float *Hg, float* x, sysstruct &sys)
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
        
//         __builtin_prefetch(&Hg[bsz2*pii+4*6250] , 0,3);
        
        // Create temporary dense vector (to perform matrix-vector product of the entire row with one BLAS call only):
        for (Int i = 0; i<nri; i++) {
            pij += 1;
            rj = col_ind[pij];
            for (Int j = 0; j < bsz; j++) {
//                 __builtin_prefetch(&x[bsz*col_ind[pij+4*250]], 0,3);
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

void matvec3(float* v, float *Hg, float* x, sysstruct &sys) {
//     DESCRIPTION: v = Hg*x, where Hg is a sparse matrix and (x,v) are vectors. Single precision implementation.

//     INPUTS / OUTPUTS:
//     Hg:      Pointer to the sparse matrix Hg.
//              Hg is stored in Block Compressed Row Storage (BCRS) format. More info: http://www.netlib.org/utk/people/JackDongarra/etemplates/node375.html
//     x:       Pointer to the input vector x.
//     v:       Pointer to the output vector v.
//     sys:     "sysstruct" structure with information about the linear system to be solved.

//     NOTES:
//     - In this function, we assume the block size is 5.

	float one = 1.0, zero = 0.0, alpha;
	char chn = 'N';
	int nri;

	Int *row_blk = &sys.row_blk[0];     // Description of this vector is in: http://www.netlib.org/utk/people/JackDongarra/etemplates/node375.html
	float *Hgpos = &Hg[0];
	float *vpos = &v[0];
	Int *col_ind_pos = &sys.col_ind[0];     // Description of this vector is in: http://www.netlib.org/utk/people/JackDongarra/etemplates/node375.html

	for (int ri = 0; ri < sys.numRows; ri++) {
		int nri = row_blk[ri+1] - row_blk[ri];

		float vtemp1 = 0;
		float vtemp2 = 0;
		float vtemp3 = 0;
		float vtemp4 = 0;
		float vtemp5 = 0;

		for (int i=0; i<nri; i++) {
			float *xpos = &x[5 * (*col_ind_pos++)];

			vtemp1 += Hgpos[0] * xpos[0] + Hgpos[5] * xpos[1] + Hgpos[10] * xpos[2] + Hgpos[15] * xpos[3] + Hgpos[20] * xpos[4];
			vtemp2 += Hgpos[1] * xpos[0] + Hgpos[6] * xpos[1] + Hgpos[11] * xpos[2] + Hgpos[16] * xpos[3] + Hgpos[21] * xpos[4];
			vtemp3 += Hgpos[2] * xpos[0] + Hgpos[7] * xpos[1] + Hgpos[12] * xpos[2] + Hgpos[17] * xpos[3] + Hgpos[22] * xpos[4];
			vtemp4 += Hgpos[3] * xpos[0] + Hgpos[8] * xpos[1] + Hgpos[13] * xpos[2] + Hgpos[18] * xpos[3] + Hgpos[23] * xpos[4];
			vtemp5 += Hgpos[4] * xpos[0] + Hgpos[9] * xpos[1] + Hgpos[14] * xpos[2] + Hgpos[19] * xpos[3] + Hgpos[24] * xpos[4];
			Hgpos += 25;
		}
		vpos[0] = vtemp1;
		vpos[1] = vtemp2;
		vpos[2] = vtemp3;
		vpos[3] = vtemp4;
		vpos[4] = vtemp5;
		vpos += 5;
	}
}

void matvec4(float* v, float *Hg, float* x, sysstruct &sys) {
//     DESCRIPTION: v = Hg*x, where Hg is a sparse matrix and (x,v) are vectors. Single precision implementation.

//     INPUTS / OUTPUTS:
//     Hg:      Pointer to the sparse matrix Hg.
//              Hg is stored in Block Compressed Row Storage (BCRS) format. More info: http://www.netlib.org/utk/people/JackDongarra/etemplates/node375.html
//     x:       Pointer to the input vector x.
//     v:       Pointer to the output vector v.
//     sys:     "sysstruct" structure with information about the linear system to be solved.

//     NOTES:
//     - In this function, we assume the block size is 5.
	
	float one = 1.0, zero = 0.0, alpha;
	char chn = 'N';
	int nri;
	
	Int *row_blk = &sys.row_blk[0];     // Description of this vector is in: http://www.netlib.org/utk/people/JackDongarra/etemplates/node375.html
	float *Hgpos = &Hg[0];
	float *vpos = &v[0];
	Int *col_ind_pos = &sys.col_ind[0];     // Description of this vector is in: http://www.netlib.org/utk/people/JackDongarra/etemplates/node375.html

	for (int ri = 0; ri < sys.numRows; ri++) {
		int nri = row_blk[ri+1] - row_blk[ri];
        
		float vtemp1 = 0;
		float vtemp2 = 0;
		float vtemp3 = 0;
		float vtemp4 = 0;
		float vtemp5 = 0;
		
		for (int i=0; i<nri; i++) {
			__builtin_prefetch(&x[5 * (*col_ind_pos + 4*250)    ] , 0,3);
			__builtin_prefetch(Hgpos+25*4*250      , 0,3);
// 			__builtin_prefetch(Hgpos+25*4*250 + 16 , 0,3);
            
			float *xpos = &x[5 * (*col_ind_pos++)];

			vtemp1 += Hgpos[0] * xpos[0] + Hgpos[5] * xpos[1] + Hgpos[10] * xpos[2] + Hgpos[15] * xpos[3] + Hgpos[20] * xpos[4];
			vtemp2 += Hgpos[1] * xpos[0] + Hgpos[6] * xpos[1] + Hgpos[11] * xpos[2] + Hgpos[16] * xpos[3] + Hgpos[21] * xpos[4];
			vtemp3 += Hgpos[2] * xpos[0] + Hgpos[7] * xpos[1] + Hgpos[12] * xpos[2] + Hgpos[17] * xpos[3] + Hgpos[22] * xpos[4];
			vtemp4 += Hgpos[3] * xpos[0] + Hgpos[8] * xpos[1] + Hgpos[13] * xpos[2] + Hgpos[18] * xpos[3] + Hgpos[23] * xpos[4];
			vtemp5 += Hgpos[4] * xpos[0] + Hgpos[9] * xpos[1] + Hgpos[14] * xpos[2] + Hgpos[19] * xpos[3] + Hgpos[24] * xpos[4];
			Hgpos += 25;
		}
		vpos[0] = vtemp1;
		vpos[1] = vtemp2;
		vpos[2] = vtemp3;
		vpos[3] = vtemp4;
		vpos[4] = vtemp5;
		vpos += 5;
	}
}
