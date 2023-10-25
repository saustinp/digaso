void fhat_ehd_tof2d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double uh1 = uh[0*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = 1.0/param2;
		double t3 = 1.0/param5;
		double t0 = nl2*x1*(u1+param1*t2*t3*u3)+param10*x1*(u1-uh1)+nl1*param1*t2*t3*u2*x1;
		fh[0*ng+i] = 0.0;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double uh1 = uh[0*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = 1.0/param2;
		double t3 = 1.0/param5;
		fh_udg[0*ng+i] = nl2*x1+param10*x1;
		fh_udg[1*ng+i] = nl1*param1*t2*t3*x1;
		fh_udg[2*ng+i] = nl2*param1*t2*t3*x1;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double uh1 = uh[0*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t0 = -param10*x1;
		fh_uh[0*ng+i] = 0.0;

	}
}

void fhatonly_ehd_tof2d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double uh1 = uh[0*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = 1.0/param2;
		double t3 = 1.0/param5;
		double t0 = nl2*x1*(u1+param1*t2*t3*u3)+param10*x1*(u1-uh1)+nl1*param1*t2*t3*u2*x1;

	}
}

