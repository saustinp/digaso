void source_streamer2d(double *s, double *s_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double t4 = t2+t3;
		double t5 = sqrt(t4);
		double t6 = 1.0/param2;
		double t7 = param3*t5;
		double t8 = 1.0/pow(t7,1.3E1/5.0E1);
		double t9 = param1*3.4075E2;
		double t10 = 1.0/param3;
		double t11 = 1.0/sqrt(t4);
		double t12 = exp(t10*t11*-2.73E7);
		double t13 = 1.0/(param3*param3*param3);
		double t14 = 1.0/pow(t4,3.0/2.0);
		double t15 = t13*t14*4.3666E26;
		double t16 = t15+1.1944E6;
		double t17 = t9-param1*t12*t16;
		s[0*ng+i] = t5*t6*t8*t17*u1*x1*(-2.3987);
		s[1*ng+i] = t5*t6*t8*t17*u1*x1*(-2.3987);
		s[2*ng+i] = -1.0/(param1*param1)*param4*param10*t10*x1*(u1-u2);

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
		double t4 = t2+t3;
		double t5 = sqrt(t4);
		double t6 = 1.0/param2;
		double t7 = param3*t5;
		double t8 = 1.0/pow(t7,1.3E1/5.0E1);
		double t9 = param1*3.4075E2;
		double t10 = 1.0/param3;
		double t11 = 1.0/sqrt(t4);
		double t19 = t10*t11*2.73E7;
		double t12 = exp(-t19);
		double t13 = 1.0/(param3*param3*param3);
		double t14 = 1.0/pow(t4,3.0/2.0);
		double t15 = t13*t14*4.3666E26;
		double t16 = t15+1.1944E6;
		double t20 = param1*t12*t16;
		double t17 = t9-t20;
		double t18 = 1.0/(param1*param1);
		double t21 = 1.0/pow(t4,5.0/2.0);
		double t22 = param1*t12*t13*t21*u6*1.30998E27;
		double t23 = t22-param1*t10*t12*t14*t16*u6*2.73E7;
		double t24 = 1.0/pow(t7,6.3E1/5.0E1);
		double t25 = param3*t6*t17*t24*u1*u6*x1*6.23662E-1;
		double t26 = t25-t5*t6*t8*t23*u1*x1*2.3987-t6*t8*t11*t17*u1*u6*x1*2.3987;
		double t27 = param1*t12*t13*t21*u9*1.30998E27;
		double t28 = t27-param1*t10*t12*t14*t16*u9*2.73E7;
		double t29 = param3*t6*t17*t24*u1*u9*x1*6.23662E-1;
		double t30 = t29-t5*t6*t8*t28*u1*x1*2.3987-t6*t8*t11*t17*u1*u9*x1*2.3987;
		s_udg[0*ng+i] = t5*t6*t8*t17*x1*(-2.3987);
		s_udg[1*ng+i] = t5*t6*t8*t17*x1*(-2.3987);
		s_udg[2*ng+i] = -param4*param10*t10*t18*x1;
		s_udg[3*ng+i] = 0.0;
		s_udg[4*ng+i] = 0.0;
		s_udg[5*ng+i] = param4*param10*t10*t18*x1;
		s_udg[6*ng+i] = 0.0;
		s_udg[7*ng+i] = 0.0;
		s_udg[8*ng+i] = 0.0;
		s_udg[9*ng+i] = 0.0;
		s_udg[10*ng+i] = 0.0;
		s_udg[11*ng+i] = 0.0;
		s_udg[12*ng+i] = 0.0;
		s_udg[13*ng+i] = 0.0;
		s_udg[14*ng+i] = 0.0;
		s_udg[15*ng+i] = t26;
		s_udg[16*ng+i] = t26;
		s_udg[17*ng+i] = 0.0;
		s_udg[18*ng+i] = 0.0;
		s_udg[19*ng+i] = 0.0;
		s_udg[20*ng+i] = 0.0;
		s_udg[21*ng+i] = 0.0;
		s_udg[22*ng+i] = 0.0;
		s_udg[23*ng+i] = 0.0;
		s_udg[24*ng+i] = t30;
		s_udg[25*ng+i] = t30;
		s_udg[26*ng+i] = 0.0;

	}
}

void sourceonly_streamer2d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double t4 = t2+t3;
		double t5 = sqrt(t4);
		double t6 = 1.0/param2;
		double t7 = param3*t5;
		double t8 = 1.0/pow(t7,1.3E1/5.0E1);
		double t9 = param1*3.4075E2;
		double t10 = 1.0/param3;
		double t11 = 1.0/sqrt(t4);
		double t12 = exp(t10*t11*-2.73E7);
		double t13 = 1.0/(param3*param3*param3);
		double t14 = 1.0/pow(t4,3.0/2.0);
		double t15 = t13*t14*4.3666E26;
		double t16 = t15+1.1944E6;
		double t17 = t9-param1*t12*t16;
		s[0*ng+i] = t5*t6*t8*t17*u1*x1*(-2.3987);
		s[1*ng+i] = t5*t6*t8*t17*u1*x1*(-2.3987);
		s[2*ng+i] = -1.0/(param1*param1)*param4*param10*t10*x1*(u1-u2);

	}
}

