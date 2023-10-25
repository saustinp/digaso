void fluxonly_ehd_tof3d(double *f, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double param4 = param[3];
	double param5 = param[4];
	double param6 = param[5];
	double param7 = param[6];
	double param8 = param[7];
	double param9 = param[8];
	double param10 = param[9];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];

		double t2 = 1.0/param2;
		double t3 = 1.0/param5;
		f[0*ng+i] = param1*t2*t3*u2*x1;
		f[1*ng+i] = x1*(u1+param1*t2*t3*u3);

	}
}

