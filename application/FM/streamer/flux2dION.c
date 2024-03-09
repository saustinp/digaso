void flux_streamer2d(double *f, double *f_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
	double param11 = param[10];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];

		double t2 = u6*u6;
		double t3 = u9*u9;
		double t4 = 1.0/param1;
		double t5 = 1.0/param2;
		double t6 = 1.0/param3;
		double t7 = t2+t3;
		double t8 = sqrt(t7);
		double t9 = param3*t8;
		double t10 = pow(t9,1.1E+1/5.0E+1);
		double t11 = 1.0/pow(t9,1.3E+1/5.0E+1);
		f[0*ng+i] = -x1*(t5*t11*u1*u6*2.3987-t4*t5*t6*t10*u4*4.3628E-3);
		f[1*ng+i] = x1*((t5*u2*u6)/5.0E+3+(t4*t5*t6*u5)/2.0E+4);
		f[2*ng+i] = u6*x1;
		f[3*ng+i] = -x1*(t5*t11*u1*u9*2.3987-t4*t5*t6*t10*u7*4.3628E-3);
		f[4*ng+i] = x1*((t5*u2*u9)/5.0E+3+(t4*t5*t6*u8)/2.0E+4);
		f[5*ng+i] = u9*x1;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];

		double t2 = u6*u6;
		double t3 = u9*u9;
		double t4 = 1.0/param1;
		double t5 = 1.0/param2;
		double t6 = 1.0/param3;
		double t7 = t2+t3;
		double t10 = (t5*u2*x1)/5.0E+3;
		double t16 = (t4*t5*t6*x1)/2.0E+4;
		double t8 = sqrt(t7);
		double t9 = 1.0/t8;
		double t11 = param3*t8;
		double t12 = pow(t11,1.1E+1/5.0E+1);
		double t13 = 1.0/pow(t11,1.3E+1/5.0E+1);
		double t15 = 1.0/pow(t11,6.3E+1/5.0E+1);
		double t14 = t13*t13*t13;
		double t17 = t5*t13*u1*2.3987;
		double t19 = param3*t5*t9*t15*u1*u6*u9*6.23662E-1;
		double t20 = t4*t5*t6*t12*x1*4.3628E-3;
		double t18 = -t17;
		f_udg[0*ng+i] = t5*t13*u6*x1*(-2.3987);
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = t5*t13*u9*x1*(-2.3987);
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = (t5*u6*x1)/5.0E+3;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = (t5*u9*x1)/5.0E+3;
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = 0.0;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = t20;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = t16;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = x1*(t18+param3*t2*t5*t9*t15*u1*6.23662E-1+t4*t5*t9*t14*u4*u6*9.59816E-4);
		f_udg[31*ng+i] = t10;
		f_udg[32*ng+i] = x1;
		f_udg[33*ng+i] = x1*(t19+t4*t5*t9*t14*u6*u7*9.59816E-4);
		f_udg[34*ng+i] = 0.0;
		f_udg[35*ng+i] = 0.0;
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = t20;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = t16;
		f_udg[47*ng+i] = 0.0;
		f_udg[48*ng+i] = x1*(t19+t4*t5*t9*t14*u4*u9*9.59816E-4);
		f_udg[49*ng+i] = 0.0;
		f_udg[50*ng+i] = 0.0;
		f_udg[51*ng+i] = x1*(t18+param3*t3*t5*t9*t15*u1*6.23662E-1+t4*t5*t9*t14*u7*u9*9.59816E-4);
		f_udg[52*ng+i] = t10;
		f_udg[53*ng+i] = x1;

	}
}

void fluxonly_streamer2d(double *f, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
	double param11 = param[10];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];

		double t2 = u6*u6;
		double t3 = u9*u9;
		double t4 = 1.0/param1;
		double t5 = 1.0/param2;
		double t6 = 1.0/param3;
		double t7 = t2+t3;
		double t8 = sqrt(t7);
		double t9 = param3*t8;
		double t10 = pow(t9,1.1E+1/5.0E+1);
		double t11 = 1.0/pow(t9,1.3E+1/5.0E+1);
		f[0*ng+i] = -x1*(t5*t11*u1*u6*2.3987-t4*t5*t6*t10*u4*4.3628E-3);
		f[1*ng+i] = x1*((t5*u2*u6)/5.0E+3+(t4*t5*t6*u5)/2.0E+4);
		f[2*ng+i] = u6*x1;
		f[3*ng+i] = -x1*(t5*t11*u1*u9*2.3987-t4*t5*t6*t10*u7*4.3628E-3);
		f[4*ng+i] = x1*((t5*u2*u9)/5.0E+3+(t4*t5*t6*u8)/2.0E+4);
		f[5*ng+i] = u9*x1;

	}
}
