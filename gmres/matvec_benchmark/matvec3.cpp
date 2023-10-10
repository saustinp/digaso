#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

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
            
            vtemp1 += (*Hgpos++) * (*xpos);
            vtemp2 += (*Hgpos++) * (*xpos);
            vtemp3 += (*Hgpos++) * (*xpos);
            vtemp4 += (*Hgpos++) * (*xpos);
            vtemp5 += (*Hgpos++) * (*xpos);

            xpos++;
            vtemp1 += (*Hgpos++) * (*xpos);
            vtemp2 += (*Hgpos++) * (*xpos);
            vtemp3 += (*Hgpos++) * (*xpos);
            vtemp4 += (*Hgpos++) * (*xpos);
            vtemp5 += (*Hgpos++) * (*xpos);

            xpos++;
            vtemp1 += (*Hgpos++) * (*xpos);
            vtemp2 += (*Hgpos++) * (*xpos);
            vtemp3 += (*Hgpos++) * (*xpos);
            vtemp4 += (*Hgpos++) * (*xpos);
            vtemp5 += (*Hgpos++) * (*xpos);

            xpos++;
            vtemp1 += (*Hgpos++) * (*xpos);
            vtemp2 += (*Hgpos++) * (*xpos);
            vtemp3 += (*Hgpos++) * (*xpos);
            vtemp4 += (*Hgpos++) * (*xpos);
            vtemp5 += (*Hgpos++) * (*xpos);

            xpos++;
            vtemp1 += (*Hgpos++) * (*xpos);
            vtemp2 += (*Hgpos++) * (*xpos);
            vtemp3 += (*Hgpos++) * (*xpos);
            vtemp4 += (*Hgpos++) * (*xpos);
            vtemp5 += (*Hgpos++) * (*xpos);
		}
		vpos[0] = vtemp1;
		vpos[1] = vtemp2;
		vpos[2] = vtemp3;
		vpos[3] = vtemp4;
		vpos[4] = vtemp5;
		vpos += 5;
	}
}
