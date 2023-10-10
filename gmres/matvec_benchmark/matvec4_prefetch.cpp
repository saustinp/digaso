#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

void matvec4_prefetch(float* v, float *Hg, float* x, sysstruct &sys, int preftchDist) {
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
			__builtin_prefetch(&x[5 * (*col_ind_pos + preftchDist)    ] , 0 , 3);
			__builtin_prefetch(Hgpos+25*preftchDist      , 0 , 2);
			__builtin_prefetch(Hgpos+25*preftchDist + 16 , 0 , 2);
            
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
