void fhat_leuq3d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double mu = param[0];
	double lambda = param[1];
	double tau = param[2];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double q11 = udg[3*ng+i];
		double q21 = udg[4*ng+i];
		double q31 = udg[5*ng+i];
		double q12 = udg[6*ng+i];
		double q22 = udg[7*ng+i];
		double q32 = udg[8*ng+i];
		double q13 = udg[9*ng+i];
		double q23 = udg[10*ng+i];
		double q33 = udg[11*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		double t2 = q11+q22+q33+3.0;
		double t3 = lambda*t2;
		double t4 = q12+q21;
		double t5 = q13+q31;
		double t6 = q23+q32;
		fh[0*ng+i] = nl1*(t3+mu*(q11+1.0)*2.0)+tau*(u1-uh1)+mu*nl2*t4+mu*nl3*t5;
		fh[1*ng+i] = nl2*(t3+mu*(q22+1.0)*2.0)+tau*(u2-uh2)+mu*nl1*t4+mu*nl3*t6;
		fh[2*ng+i] = nl3*(t3+mu*(q33+1.0)*2.0)+tau*(u3-uh3)+mu*nl1*t5+mu*nl2*t6;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double q11 = udg[3*ng+i];
		double q21 = udg[4*ng+i];
		double q31 = udg[5*ng+i];
		double q12 = udg[6*ng+i];
		double q22 = udg[7*ng+i];
		double q32 = udg[8*ng+i];
		double q13 = udg[9*ng+i];
		double q23 = udg[10*ng+i];
		double q33 = udg[11*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		double t2 = mu*nl1;
		double t3 = mu*nl2;
		double t4 = mu*2.0;
		double t5 = lambda+t4;
		double t6 = lambda*nl3;
		double t7 = mu*nl3;
		double t8 = lambda*nl1;
		double t9 = lambda*nl2;
		fh_udg[0*ng+i] = tau;
		fh_udg[1*ng+i] = 0.0;
		fh_udg[2*ng+i] = 0.0;
		fh_udg[3*ng+i] = 0.0;
		fh_udg[4*ng+i] = tau;
		fh_udg[5*ng+i] = 0.0;
		fh_udg[6*ng+i] = 0.0;
		fh_udg[7*ng+i] = 0.0;
		fh_udg[8*ng+i] = tau;
		fh_udg[9*ng+i] = nl1*t5;
		fh_udg[10*ng+i] = t9;
		fh_udg[11*ng+i] = t6;
		fh_udg[12*ng+i] = t3;
		fh_udg[13*ng+i] = t2;
		fh_udg[14*ng+i] = 0.0;
		fh_udg[15*ng+i] = t7;
		fh_udg[16*ng+i] = 0.0;
		fh_udg[17*ng+i] = t2;
		fh_udg[18*ng+i] = t3;
		fh_udg[19*ng+i] = t2;
		fh_udg[20*ng+i] = 0.0;
		fh_udg[21*ng+i] = t8;
		fh_udg[22*ng+i] = nl2*t5;
		fh_udg[23*ng+i] = t6;
		fh_udg[24*ng+i] = 0.0;
		fh_udg[25*ng+i] = t7;
		fh_udg[26*ng+i] = t3;
		fh_udg[27*ng+i] = t7;
		fh_udg[28*ng+i] = 0.0;
		fh_udg[29*ng+i] = t2;
		fh_udg[30*ng+i] = 0.0;
		fh_udg[31*ng+i] = t7;
		fh_udg[32*ng+i] = t3;
		fh_udg[33*ng+i] = t8;
		fh_udg[34*ng+i] = t9;
		fh_udg[35*ng+i] = nl3*t5;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double q11 = udg[3*ng+i];
		double q21 = udg[4*ng+i];
		double q31 = udg[5*ng+i];
		double q12 = udg[6*ng+i];
		double q22 = udg[7*ng+i];
		double q32 = udg[8*ng+i];
		double q13 = udg[9*ng+i];
		double q23 = udg[10*ng+i];
		double q33 = udg[11*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		fh_uh[0*ng+i] = -tau;
		fh_uh[1*ng+i] = 0.0;
		fh_uh[2*ng+i] = 0.0;
		fh_uh[3*ng+i] = 0.0;
		fh_uh[4*ng+i] = -tau;
		fh_uh[5*ng+i] = 0.0;
		fh_uh[6*ng+i] = 0.0;
		fh_uh[7*ng+i] = 0.0;
		fh_uh[8*ng+i] = -tau;

	}
}

void fhatonly_leuq3d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double mu = param[0];
	double lambda = param[1];
	double tau = param[2];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double q11 = udg[3*ng+i];
		double q21 = udg[4*ng+i];
		double q31 = udg[5*ng+i];
		double q12 = udg[6*ng+i];
		double q22 = udg[7*ng+i];
		double q32 = udg[8*ng+i];
		double q13 = udg[9*ng+i];
		double q23 = udg[10*ng+i];
		double q33 = udg[11*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		double t2 = q11+q22+q33+3.0;
		double t3 = lambda*t2;
		double t4 = q12+q21;
		double t5 = q13+q31;
		double t6 = q23+q32;
		fh[0*ng+i] = nl1*(t3+mu*(q11+1.0)*2.0)+tau*(u1-uh1)+mu*nl2*t4+mu*nl3*t5;
		fh[1*ng+i] = nl2*(t3+mu*(q22+1.0)*2.0)+tau*(u2-uh2)+mu*nl1*t4+mu*nl3*t6;
		fh[2*ng+i] = nl3*(t3+mu*(q33+1.0)*2.0)+tau*(u3-uh3)+mu*nl1*t5+mu*nl2*t6;

	}
}

