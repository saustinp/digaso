void sourceonly_axispoisson3d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	int i;
	for (i = 0; i <ng*ncu; i++)
		s[i] = 0.0;
}

