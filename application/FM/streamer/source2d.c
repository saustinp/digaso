// THIS IS THE LOGARITHMIC FORMULATION!

void source2d(double *s, double *s_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];

		double t2 = exp(u1);
		double t3 = exp(u2);
		double t4 = u6*u6;
		double t5 = u9*u9;
		double t6 = 1.0/(param1*param1);
		double t7 = 1.0/param2;
		double t8 = 1.0/param3;
		double t13 = param1*3.4075E+2;
		double t9 = t8*t8*t8;
		double t10 = -t3;
		double t11 = t4+t5;
		double t12 = t2+t10;
		double t14 = sqrt(t11);
		double t15 = 1.0/t14;
		double t17 = param3*t14;
		double t16 = t15*t15*t15;
		double t18 = 1.0/pow(t17,1.3E+1/5.0E+1);
		double t19 = t8*t15*2.73E+7;
		double t20 = -t19;
		double t22 = t9*t16*4.3666E+26;
		double t21 = exp(t20);
		double t23 = t22+1.1944E+6;
		double t24 = param1*t21*t23;
		double t25 = -t24;
		double t26 = t13+t25;
		s[0*ng+i] = x1*(t7*t14*t18*t26*(-2.3987)+(t7*t8*pow(t17,1.1E+1/5.0E+1)*(u4*u4+u7*u7)*4.3628E-3)/param1+param4*t6*t7*t8*t12*t18*(u1-1.0)*2.3987);
		s[1*ng+i] = t7*t14*t18*t26*x1*exp(u1-u2)*(-2.3987);
		s[2*ng+i] = -param4*t6*t8*t12*x1;

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

		double t2 = exp(u1);
		double t3 = exp(u2);
		double t4 = u4*u4;
		double t5 = u6*u6;
		double t6 = u7*u7;
		double t7 = u9*u9;
		double t8 = 1.0/param1;
		double t10 = 1.0/param2;
		double t11 = 1.0/param3;
		double t13 = -u2;
		double t14 = u1-1.0;
		double t21 = param1*3.4075E+2;
		double t9 = t8*t8;
		double t12 = t11*t11*t11;
		double t15 = -t3;
		double t16 = t13+u1;
		double t17 = t4+t6;
		double t18 = t5+t7;
		double t19 = exp(t16);
		double t20 = t2+t15;
		double t22 = sqrt(t18);
		double t23 = 1.0/t22;
		double t26 = param3*t22;
		double t24 = t23*t23*t23;
		double t25 = t23*t23*t23*t23*t23;
		double t27 = pow(t26,1.1E+1/5.0E+1);
		double t28 = 1.0/pow(t26,1.3E+1/5.0E+1);
		double t30 = 1.0/pow(t26,6.3E+1/5.0E+1);
		double t31 = t11*t23*2.73E+7;
		double t29 = t28*t28*t28;
		double t32 = -t31;
		double t34 = t12*t24*4.3666E+26;
		double t33 = exp(t32);
		double t35 = t34+1.1944E+6;
		double t36 = param1*t12*t25*t33*u6*1.30998E+27;
		double t37 = param1*t12*t25*t33*u9*1.30998E+27;
		double t38 = param1*t33*t35;
		double t39 = -t38;
		double t41 = t11*t24*t38*u6*2.73E+7;
		double t42 = t11*t24*t38*u9*2.73E+7;
		double t40 = t21+t39;
		double t43 = -t41;
		double t44 = -t42;
		double t45 = t10*t19*t22*t28*t40*x1*2.3987;
		double t46 = t36+t43;
		double t47 = t37+t44;
		s_udg[0*ng+i] = x1*(param4*t9*t10*t11*t20*t28*2.3987+param4*t2*t9*t10*t11*t14*t28*2.3987);
		s_udg[1*ng+i] = -t45;
		s_udg[2*ng+i] = -param4*t2*t9*t11*x1;
		s_udg[3*ng+i] = param4*t3*t9*t10*t11*t14*t28*x1*(-2.3987);
		s_udg[4*ng+i] = t45;
		s_udg[5*ng+i] = param4*t3*t9*t11*x1;
		s_udg[6*ng+i] = 0.0;
		s_udg[7*ng+i] = 0.0;
		s_udg[8*ng+i] = 0.0;
		s_udg[9*ng+i] = t8*t10*t11*t27*u4*x1*8.7256E-3;
		s_udg[10*ng+i] = 0.0;
		s_udg[11*ng+i] = 0.0;
		s_udg[12*ng+i] = 0.0;
		s_udg[13*ng+i] = 0.0;
		s_udg[14*ng+i] = 0.0;
		s_udg[15*ng+i] = -x1*(t10*t22*t28*t46*2.3987-param3*t10*t30*t40*u6*6.23662E-1+t10*t23*t28*t40*u6*2.3987-t8*t10*t17*t23*t29*u6*9.59816E-4+param4*t9*t10*t14*t20*t23*t30*u6*6.23662E-1);
		s_udg[16*ng+i] = t10*t19*t22*t28*t46*x1*(-2.3987)+param3*t10*t19*t30*t40*u6*x1*6.23662E-1-t10*t19*t23*t28*t40*u6*x1*2.3987;
		s_udg[17*ng+i] = 0.0;
		s_udg[18*ng+i] = t8*t10*t11*t27*u7*x1*8.7256E-3;
		s_udg[19*ng+i] = 0.0;
		s_udg[20*ng+i] = 0.0;
		s_udg[21*ng+i] = 0.0;
		s_udg[22*ng+i] = 0.0;
		s_udg[23*ng+i] = 0.0;
		s_udg[24*ng+i] = -x1*(t10*t22*t28*t47*2.3987-param3*t10*t30*t40*u9*6.23662E-1+t10*t23*t28*t40*u9*2.3987-t8*t10*t17*t23*t29*u9*9.59816E-4+param4*t9*t10*t14*t20*t23*t30*u9*6.23662E-1);
		s_udg[25*ng+i] = t10*t19*t22*t28*t47*x1*(-2.3987)+param3*t10*t19*t30*t40*u9*x1*6.23662E-1-t10*t19*t23*t28*t40*u9*x1*2.3987;
		s_udg[26*ng+i] = 0.0;

	}
}

void sourceonly2d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];

		double t2 = exp(u1);
		double t3 = exp(u2);
		double t4 = u6*u6;
		double t5 = u9*u9;
		double t6 = 1.0/(param1*param1);
		double t7 = 1.0/param2;
		double t8 = 1.0/param3;
		double t13 = param1*3.4075E+2;
		double t9 = t8*t8*t8;
		double t10 = -t3;
		double t11 = t4+t5;
		double t12 = t2+t10;
		double t14 = sqrt(t11);
		double t15 = 1.0/t14;
		double t17 = param3*t14;
		double t16 = t15*t15*t15;
		double t18 = 1.0/pow(t17,1.3E+1/5.0E+1);
		double t19 = t8*t15*2.73E+7;
		double t20 = -t19;
		double t22 = t9*t16*4.3666E+26;
		double t21 = exp(t20);
		double t23 = t22+1.1944E+6;
		double t24 = param1*t21*t23;
		double t25 = -t24;
		double t26 = t13+t25;
		s[0*ng+i] = x1*(t7*t14*t18*t26*(-2.3987)+(t7*t8*pow(t17,1.1E+1/5.0E+1)*(u4*u4+u7*u7)*4.3628E-3)/param1+param4*t6*t7*t8*t12*t18*(u1-1.0)*2.3987);
		s[1*ng+i] = t7*t14*t18*t26*x1*exp(u1-u2)*(-2.3987);
		s[2*ng+i] = -param4*t6*t8*t12*x1;

	}
}

