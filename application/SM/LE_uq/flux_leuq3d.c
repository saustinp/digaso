void flux_leuq3d(double *f, double *f_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];

    //cout<<"flux"<<endl;
    
	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];
		double u10 = udg[9*ng+i];
		double u11 = udg[10*ng+i];
		double u12 = udg[11*ng+i];

		double t2 = u5+u7;
		double t3 = param1*t2;
		double t4 = u4+u8+u12+3.0;
		double t5 = param2*t4;
		double t6 = u6+u10;
		double t7 = param1*t6;
		double t8 = u9+u11;
		double t9 = param1*t8;
		f[0*ng+i] = t5+param1*(u4+1.0)*2.0;
		f[1*ng+i] = t3;
		f[2*ng+i] = t7;
		f[3*ng+i] = t3;
		f[4*ng+i] = t5+param1*(u8+1.0)*2.0;
		f[5*ng+i] = t9;
		f[6*ng+i] = t7;
		f[7*ng+i] = t9;
		f[8*ng+i] = t5+param1*(u12+1.0)*2.0;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];
		double u10 = udg[9*ng+i];
		double u11 = udg[10*ng+i];
		double u12 = udg[11*ng+i];

		double t2 = param1*2.0;
		double t3 = param2+t2;
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = 0.0;
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = 0.0;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = 0.0;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = 0.0;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = t3;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = param2;
		f_udg[32*ng+i] = 0.0;
		f_udg[33*ng+i] = 0.0;
		f_udg[34*ng+i] = 0.0;
		f_udg[35*ng+i] = param2;
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = param1;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = param1;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = param1;
		f_udg[48*ng+i] = 0.0;
		f_udg[49*ng+i] = 0.0;
		f_udg[50*ng+i] = 0.0;
		f_udg[51*ng+i] = param1;
		f_udg[52*ng+i] = 0.0;
		f_udg[53*ng+i] = 0.0;
		f_udg[54*ng+i] = 0.0;
		f_udg[55*ng+i] = param1;
		f_udg[56*ng+i] = 0.0;
		f_udg[57*ng+i] = param1;
		f_udg[58*ng+i] = 0.0;
		f_udg[59*ng+i] = 0.0;
		f_udg[60*ng+i] = 0.0;
		f_udg[61*ng+i] = 0.0;
		f_udg[62*ng+i] = 0.0;
		f_udg[63*ng+i] = param2;
		f_udg[64*ng+i] = 0.0;
		f_udg[65*ng+i] = 0.0;
		f_udg[66*ng+i] = 0.0;
		f_udg[67*ng+i] = t3;
		f_udg[68*ng+i] = 0.0;
		f_udg[69*ng+i] = 0.0;
		f_udg[70*ng+i] = 0.0;
		f_udg[71*ng+i] = param2;
		f_udg[72*ng+i] = 0.0;
		f_udg[73*ng+i] = 0.0;
		f_udg[74*ng+i] = 0.0;
		f_udg[75*ng+i] = 0.0;
		f_udg[76*ng+i] = 0.0;
		f_udg[77*ng+i] = param1;
		f_udg[78*ng+i] = 0.0;
		f_udg[79*ng+i] = param1;
		f_udg[80*ng+i] = 0.0;
		f_udg[81*ng+i] = 0.0;
		f_udg[82*ng+i] = 0.0;
		f_udg[83*ng+i] = param1;
		f_udg[84*ng+i] = 0.0;
		f_udg[85*ng+i] = 0.0;
		f_udg[86*ng+i] = 0.0;
		f_udg[87*ng+i] = param1;
		f_udg[88*ng+i] = 0.0;
		f_udg[89*ng+i] = 0.0;
		f_udg[90*ng+i] = 0.0;
		f_udg[91*ng+i] = 0.0;
		f_udg[92*ng+i] = 0.0;
		f_udg[93*ng+i] = 0.0;
		f_udg[94*ng+i] = 0.0;
		f_udg[95*ng+i] = param1;
		f_udg[96*ng+i] = 0.0;
		f_udg[97*ng+i] = param1;
		f_udg[98*ng+i] = 0.0;
		f_udg[99*ng+i] = param2;
		f_udg[100*ng+i] = 0.0;
		f_udg[101*ng+i] = 0.0;
		f_udg[102*ng+i] = 0.0;
		f_udg[103*ng+i] = param2;
		f_udg[104*ng+i] = 0.0;
		f_udg[105*ng+i] = 0.0;
		f_udg[106*ng+i] = 0.0;
		f_udg[107*ng+i] = t3;

	}
}

void fluxonly_leuq3d(double *f, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];
		double u10 = udg[9*ng+i];
		double u11 = udg[10*ng+i];
		double u12 = udg[11*ng+i];

		double t2 = u5+u7;
		double t3 = param1*t2;
		double t4 = u4+u8+u12+3.0;
		double t5 = param2*t4;
		double t6 = u6+u10;
		double t7 = param1*t6;
		double t8 = u9+u11;
		double t9 = param1*t8;
		f[0*ng+i] = t5+param1*(u4+1.0)*2.0;
		f[1*ng+i] = t3;
		f[2*ng+i] = t7;
		f[3*ng+i] = t3;
		f[4*ng+i] = t5+param1*(u8+1.0)*2.0;
		f[5*ng+i] = t9;
		f[6*ng+i] = t7;
		f[7*ng+i] = t9;
		f[8*ng+i] = t5+param1*(u12+1.0)*2.0;

	}
}

