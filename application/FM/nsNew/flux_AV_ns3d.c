
// Written by C. Nguyen and P. Fernandez

void flux_AV_ns3d(double *f, double *f_udg, double *pg, double *udg, appstruct &app, double *param, double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
// Artificial viscosity fluxes. Latest version of the model.
		
	double rampFactor = app.rampFactor;
	double porder = (double) app.porder[0];
	double alpha = 100.0;
	double beta = 1.0e-2;
	double k_h = 1.5;
	double eps_v = 1.0e-8;

	double gam = param[0];
	double Minf = param[4];

	for (int i = 0; i <ng; i++) {
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
		double u13 = udg[12*ng+i];
		double u14 = udg[13*ng+i];
		double u15 = udg[14*ng+i];
		double u16 = udg[15*ng+i];
		double u17 = udg[16*ng+i];
		double u18 = udg[17*ng+i];
		double u19 = udg[18*ng+i];
		double u20 = udg[19*ng+i];

		double h = pg[3*ng+i];

		double r = u1;
		double ru = u2;
		double rv = u3;
		double rw = u4;
		double rE = u5;
		double rx = u6;
		double rux = u7;
		double rvx = u8;
		double rwx = u9;
		double rEx = u10;
		double ry = u11;
		double ruy = u12;
		double rvy = u13;
		double rwy = u14;
		double rEy = u15;
		double rz = u16;
		double ruz = u17;
		double rvz = u18;
		double rwz = u19;
		double rEz = u20;

		double u = ru/r;
		double v = rv/r;
		double w = rw/r;
		double ux = (rux - u*rx)/r;
		double vy = (rvy - v*ry)/r;
		double wz = (rwz - w*rz)/r;

		double t2 = 1.0/porder;
		double t3 = 1.0/r;
		double t15 = ru*rx*t3;
		double t4 = rux-t15;
		double t5 = t3*t4;
		double t16 = rv*ry*t3;
		double t6 = rvy-t16;
		double t7 = t3*t6;
		double t17 = rw*rz*t3;
		double t8 = rwz-t17;
		double t9 = t3*t8;
		double t10 = t5+t7+t9;
		double t18 = h*k_h*t2*t10;
		double t11 = beta-t18;
		double t12 = 1.0/3.141592653589793;
		double t13 = 1.0/Minf;
		double t14 = t13+1.0;
		double t19 = alpha*t11;
		double t20 = atan(t19);
		double t21 = t12*t20;
		double t22 = t21-1.0/2.0;
		double t23 = t11*t22;
		double t24 = atan(alpha);
		double t26 = t12*t24;
		double t25 = t23-t26+1.0/2.0;
		double t27 = 1.0/(r*r);
		double t28 = gam*(1.0/2.0);
		double t29 = t28-1.0/2.0;
		double t30 = ru*ru;
		double t31 = rv*rv;
		double t32 = rw*rw;
		f[0*ng+i] += h*k_h*rampFactor*rx*t2*t14*t25;
		f[1*ng+i] += h*k_h*rampFactor*rux*t2*t14*t25;
		f[2*ng+i] += h*k_h*rampFactor*rvx*t2*t14*t25;
		f[3*ng+i] += h*k_h*rampFactor*rwx*t2*t14*t25;
		f[4*ng+i] += h*k_h*rampFactor*t2*t14*t25*(gam*rEx-t29*(ru*t3*t4*2.0+rx*t27*t30+rx*t27*t31+rx*t27*t32+rv*t3*(rvx-rv*rx*t3)*2.0+rw*t3*(rwx-rw*rx*t3)*2.0));
		f[5*ng+i] += h*k_h*rampFactor*ry*t2*t14*t25;
		f[6*ng+i] += h*k_h*rampFactor*ruy*t2*t14*t25;
		f[7*ng+i] += h*k_h*rampFactor*rvy*t2*t14*t25;
		f[8*ng+i] += h*k_h*rampFactor*rwy*t2*t14*t25;
		f[9*ng+i] += h*k_h*rampFactor*t2*t14*t25*(gam*rEy-t29*(rv*t3*t6*2.0+ry*t27*t30+ry*t27*t31+ry*t27*t32+ru*t3*(ruy-ru*ry*t3)*2.0+rw*t3*(rwy-rw*ry*t3)*2.0));
		f[10*ng+i] += h*k_h*rampFactor*rz*t2*t14*t25;
		f[11*ng+i] += h*k_h*rampFactor*ruz*t2*t14*t25;
		f[12*ng+i] += h*k_h*rampFactor*rvz*t2*t14*t25;
		f[13*ng+i] += h*k_h*rampFactor*rwz*t2*t14*t25;
		f[14*ng+i] += h*k_h*rampFactor*t2*t14*t25*(gam*rEz-t29*(rw*t3*t8*2.0+rz*t27*t30+rz*t27*t31+rz*t27*t32+ru*t3*(ruz-ru*rz*t3)*2.0+rv*t3*(rvz-rv*rz*t3)*2.0));

	}

	if (computeJacobian == 1) {

		for (int i = 0; i <ng; i++) {
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
			double u13 = udg[12*ng+i];
			double u14 = udg[13*ng+i];
			double u15 = udg[14*ng+i];
			double u16 = udg[15*ng+i];
			double u17 = udg[16*ng+i];
			double u18 = udg[17*ng+i];
			double u19 = udg[18*ng+i];
			double u20 = udg[19*ng+i];

			double h = pg[3*ng+i];

            double r = u1;
            double ru = u2;
            double rv = u3;
            double rw = u4;
            double rE = u5;
            double rx = u6;
            double rux = u7;
            double rvx = u8;
            double rwx = u9;
            double rEx = u10;
            double ry = u11;
            double ruy = u12;
            double rvy = u13;
            double rwy = u14;
            double rEy = u15;
            double rz = u16;
            double ruz = u17;
            double rvz = u18;
            double rwz = u19;
            double rEz = u20;

            double u = ru/r;
            double v = rv/r;
            double w = rw/r;
            double ux = (rux - u*rx)/r;
            double vy = (rvy - v*ry)/r;
            double wz = (rwz - w*rz)/r;

			double t2 = 1.0/porder;
			double t3 = 1.0/r;
			double t10 = ru*rx*t3;
			double t4 = rux-t10;
			double t5 = 1.0/(r*r);
			double t12 = rv*ry*t3;
			double t6 = rvy-t12;
			double t14 = rw*rz*t3;
			double t7 = rwz-t14;
			double t8 = 1.0/(r*r*r);
			double t9 = 1.0/3.141592653589793;
			double t11 = t3*t4;
			double t13 = t3*t6;
			double t15 = t3*t7;
			double t16 = t11+t13+t15;
			double t18 = h*k_h*t2*t16;
			double t17 = beta-t18;
			double t19 = t4*t5;
			double t20 = t5*t6;
			double t21 = t5*t7;
			double t29 = ru*rx*t8;
			double t30 = rv*ry*t8;
			double t31 = rw*rz*t8;
			double t22 = t19+t20+t21-t29-t30-t31;
			double t23 = 1.0/Minf;
			double t24 = t23+1.0;
			double t25 = alpha*t17;
			double t26 = atan(t25);
			double t27 = t9*t26;
			double t28 = t27-1.0/2.0;
			double t32 = h*k_h*t2*t22*t28;
			double t33 = alpha*alpha;
			double t34 = t17*t17;
			double t35 = t33*t34;
			double t36 = t35+1.0;
			double t37 = 1.0/t36;
			double t38 = alpha*h*k_h*t2*t9*t17*t22*t37;
			double t39 = t32+t38;
			double t40 = gam*(1.0/2.0);
			double t41 = t40-1.0/2.0;
			double t63 = rv*rx*t3;
			double t42 = rvx-t63;
			double t65 = rw*rx*t3;
			double t43 = rwx-t65;
			double t44 = ru*ru;
			double t45 = rv*rv;
			double t46 = rw*rw;
			double t73 = ru*ry*t3;
			double t47 = ruy-t73;
			double t76 = rw*ry*t3;
			double t48 = rwy-t76;
			double t49 = t17*t28;
			double t50 = atan(alpha);
			double t54 = t9*t50;
			double t51 = t49-t54+1.0/2.0;
			double t84 = ru*rz*t3;
			double t52 = ruz-t84;
			double t86 = rv*rz*t3;
			double t53 = rvz-t86;
			double t55 = h*k_h*rx*t2*t5*t28;
			double t56 = alpha*h*k_h*rx*t2*t5*t9*t17*t37;
			double t57 = t55+t56;
			double t58 = gam*rEx;
			double t59 = rx*t5*t44;
			double t60 = rx*t5*t45;
			double t61 = rx*t5*t46;
			double t62 = ru*t3*t4*2.0;
			double t64 = rv*t3*t42*2.0;
			double t66 = rw*t3*t43*2.0;
			double t67 = t59+t60+t61+t62+t64+t66;
			double t94 = t41*t67;
			double t68 = t58-t94;
			double t69 = gam*rEy;
			double t70 = ry*t5*t44;
			double t71 = ry*t5*t45;
			double t72 = ry*t5*t46;
			double t74 = ru*t3*t47*2.0;
			double t75 = rv*t3*t6*2.0;
			double t77 = rw*t3*t48*2.0;
			double t78 = t70+t71+t72+t74+t75+t77;
			double t95 = t41*t78;
			double t79 = t69-t95;
			double t80 = gam*rEz;
			double t81 = rz*t5*t44;
			double t82 = rz*t5*t45;
			double t83 = rz*t5*t46;
			double t85 = ru*t3*t52*2.0;
			double t87 = rv*t3*t53*2.0;
			double t88 = rw*t3*t7*2.0;
			double t89 = t81+t82+t83+t85+t87+t88;
			double t96 = t41*t89;
			double t90 = t80-t96;
			double t91 = h*k_h*ry*t2*t5*t28;
			double t92 = alpha*h*k_h*ry*t2*t5*t9*t17*t37;
			double t93 = t91+t92;
			double t97 = h*k_h*rz*t2*t5*t28;
			double t98 = alpha*h*k_h*rz*t2*t5*t9*t17*t37;
			double t99 = t97+t98;
			double t100 = h*k_h*ru*t2*t5*t28;
			double t101 = alpha*h*k_h*ru*t2*t5*t9*t17*t37;
			double t102 = t100+t101;
			double t103 = h*k_h*rampFactor*t2*t24*t51;
			double t104 = h*k_h*t2*t3*t28;
			double t105 = alpha*h*k_h*t2*t3*t9*t17*t37;
			double t106 = t104+t105;
			double t107 = h*k_h*rv*t2*t5*t28;
			double t108 = alpha*h*k_h*rv*t2*t5*t9*t17*t37;
			double t109 = t107+t108;
			double t110 = t5*t44;
			double t111 = t5*t45;
			double t112 = t5*t46;
			double t113 = t110+t111+t112;
			double t114 = h*k_h*rampFactor*t2*t24*t41*t51*t113;
			double t115 = gam*h*k_h*rampFactor*t2*t24*t51;
			double t116 = h*k_h*rw*t2*t5*t28;
			double t117 = alpha*h*k_h*rw*t2*t5*t9*t17*t37;
			double t118 = t116+t117;
			f_udg[0*ng+i] += h*k_h*rampFactor*rx*t2*t24*t39;
			f_udg[1*ng+i] += h*k_h*rampFactor*rux*t2*t24*t39;
			f_udg[2*ng+i] += h*k_h*rampFactor*rvx*t2*t24*t39;
			f_udg[3*ng+i] += h*k_h*rampFactor*rwx*t2*t24*t39;
			f_udg[4*ng+i] += h*k_h*rampFactor*t2*t24*t39*t68+h*k_h*rampFactor*t2*t24*t41*t51*(ru*t4*t5*2.0+rv*t5*t42*2.0+rw*t5*t43*2.0);
			f_udg[5*ng+i] += h*k_h*rampFactor*ry*t2*t24*t39;
			f_udg[6*ng+i] += h*k_h*rampFactor*ruy*t2*t24*t39;
			f_udg[7*ng+i] += h*k_h*rampFactor*rvy*t2*t24*t39;
			f_udg[8*ng+i] += h*k_h*rampFactor*rwy*t2*t24*t39;
			f_udg[9*ng+i] += h*k_h*rampFactor*t2*t24*t39*t79+h*k_h*rampFactor*t2*t24*t41*t51*(ru*t5*t47*2.0+rv*t5*t6*2.0+rw*t5*t48*2.0);
			f_udg[10*ng+i] += h*k_h*rampFactor*rz*t2*t24*t39;
			f_udg[11*ng+i] += h*k_h*rampFactor*ruz*t2*t24*t39;
			f_udg[12*ng+i] += h*k_h*rampFactor*rvz*t2*t24*t39;
			f_udg[13*ng+i] += h*k_h*rampFactor*rwz*t2*t24*t39;
			f_udg[14*ng+i] += h*k_h*rampFactor*t2*t24*t39*t90+h*k_h*rampFactor*t2*t24*t41*t51*(ru*t5*t52*2.0+rv*t5*t53*2.0+rw*t5*t7*2.0);
			f_udg[15*ng+i] += h*k_h*rampFactor*rx*t2*t24*t57;
			f_udg[16*ng+i] += h*k_h*rampFactor*rux*t2*t24*t57;
			f_udg[17*ng+i] += h*k_h*rampFactor*rvx*t2*t24*t57;
			f_udg[18*ng+i] += h*k_h*rampFactor*rwx*t2*t24*t57;
			f_udg[19*ng+i] += h*k_h*rampFactor*t2*t24*t57*t68-h*k_h*rampFactor*t2*t3*t4*t24*t41*t51*2.0;
			f_udg[20*ng+i] += h*k_h*rampFactor*ry*t2*t24*t57;
			f_udg[21*ng+i] += h*k_h*rampFactor*ruy*t2*t24*t57;
			f_udg[22*ng+i] += h*k_h*rampFactor*rvy*t2*t24*t57;
			f_udg[23*ng+i] += h*k_h*rampFactor*rwy*t2*t24*t57;
			f_udg[24*ng+i] += h*k_h*rampFactor*t2*t24*t57*t79-h*k_h*rampFactor*t2*t3*t24*t41*t47*t51*2.0;
			f_udg[25*ng+i] += h*k_h*rampFactor*rz*t2*t24*t57;
			f_udg[26*ng+i] += h*k_h*rampFactor*ruz*t2*t24*t57;
			f_udg[27*ng+i] += h*k_h*rampFactor*rvz*t2*t24*t57;
			f_udg[28*ng+i] += h*k_h*rampFactor*rwz*t2*t24*t57;
			f_udg[29*ng+i] += h*k_h*rampFactor*t2*t24*t57*t90-h*k_h*rampFactor*t2*t3*t24*t41*t51*t52*2.0;
			f_udg[30*ng+i] += h*k_h*rampFactor*rx*t2*t24*t93;
			f_udg[31*ng+i] += h*k_h*rampFactor*rux*t2*t24*t93;
			f_udg[32*ng+i] += h*k_h*rampFactor*rvx*t2*t24*t93;
			f_udg[33*ng+i] += h*k_h*rampFactor*rwx*t2*t24*t93;
			f_udg[34*ng+i] += h*k_h*rampFactor*t2*t24*t68*t93-h*k_h*rampFactor*t2*t3*t24*t41*t42*t51*2.0;
			f_udg[35*ng+i] += h*k_h*rampFactor*ry*t2*t24*t93;
			f_udg[36*ng+i] += h*k_h*rampFactor*ruy*t2*t24*t93;
			f_udg[37*ng+i] += h*k_h*rampFactor*rvy*t2*t24*t93;
			f_udg[38*ng+i] += h*k_h*rampFactor*rwy*t2*t24*t93;
			f_udg[39*ng+i] += h*k_h*rampFactor*t2*t24*t79*t93-h*k_h*rampFactor*t2*t3*t6*t24*t41*t51*2.0;
			f_udg[40*ng+i] += h*k_h*rampFactor*rz*t2*t24*t93;
			f_udg[41*ng+i] += h*k_h*rampFactor*ruz*t2*t24*t93;
			f_udg[42*ng+i] += h*k_h*rampFactor*rvz*t2*t24*t93;
			f_udg[43*ng+i] += h*k_h*rampFactor*rwz*t2*t24*t93;
			f_udg[44*ng+i] += h*k_h*rampFactor*t2*t24*t90*t93-h*k_h*rampFactor*t2*t3*t24*t41*t51*t53*2.0;
			f_udg[45*ng+i] += h*k_h*rampFactor*rx*t2*t24*t99;
			f_udg[46*ng+i] += h*k_h*rampFactor*rux*t2*t24*t99;
			f_udg[47*ng+i] += h*k_h*rampFactor*rvx*t2*t24*t99;
			f_udg[48*ng+i] += h*k_h*rampFactor*rwx*t2*t24*t99;
			f_udg[49*ng+i] += h*k_h*rampFactor*t2*t24*t68*t99-h*k_h*rampFactor*t2*t3*t24*t41*t43*t51*2.0;
			f_udg[50*ng+i] += h*k_h*rampFactor*ry*t2*t24*t99;
			f_udg[51*ng+i] += h*k_h*rampFactor*ruy*t2*t24*t99;
			f_udg[52*ng+i] += h*k_h*rampFactor*rvy*t2*t24*t99;
			f_udg[53*ng+i] += h*k_h*rampFactor*rwy*t2*t24*t99;
			f_udg[54*ng+i] += h*k_h*rampFactor*t2*t24*t79*t99-h*k_h*rampFactor*t2*t3*t24*t41*t48*t51*2.0;
			f_udg[55*ng+i] += h*k_h*rampFactor*rz*t2*t24*t99;
			f_udg[56*ng+i] += h*k_h*rampFactor*ruz*t2*t24*t99;
			f_udg[57*ng+i] += h*k_h*rampFactor*rvz*t2*t24*t99;
			f_udg[58*ng+i] += h*k_h*rampFactor*rwz*t2*t24*t99;
			f_udg[59*ng+i] += h*k_h*rampFactor*t2*t24*t90*t99-h*k_h*rampFactor*t2*t3*t7*t24*t41*t51*2.0;
			f_udg[60*ng+i] += 0.0;
			f_udg[61*ng+i] += 0.0;
			f_udg[62*ng+i] += 0.0;
			f_udg[63*ng+i] += 0.0;
			f_udg[64*ng+i] += 0.0;
			f_udg[65*ng+i] += 0.0;
			f_udg[66*ng+i] += 0.0;
			f_udg[67*ng+i] += 0.0;
			f_udg[68*ng+i] += 0.0;
			f_udg[69*ng+i] += 0.0;
			f_udg[70*ng+i] += 0.0;
			f_udg[71*ng+i] += 0.0;
			f_udg[72*ng+i] += 0.0;
			f_udg[73*ng+i] += 0.0;
			f_udg[74*ng+i] += 0.0;
			f_udg[75*ng+i] += t103+h*k_h*rampFactor*rx*t2*t24*t102;
			f_udg[76*ng+i] += h*k_h*rampFactor*rux*t2*t24*t102;
			f_udg[77*ng+i] += h*k_h*rampFactor*rvx*t2*t24*t102;
			f_udg[78*ng+i] += h*k_h*rampFactor*rwx*t2*t24*t102;
			f_udg[79*ng+i] += t114+h*k_h*rampFactor*t2*t24*t68*t102;
			f_udg[80*ng+i] += h*k_h*rampFactor*ry*t2*t24*t102;
			f_udg[81*ng+i] += h*k_h*rampFactor*ruy*t2*t24*t102;
			f_udg[82*ng+i] += h*k_h*rampFactor*rvy*t2*t24*t102;
			f_udg[83*ng+i] += h*k_h*rampFactor*rwy*t2*t24*t102;
			f_udg[84*ng+i] += h*k_h*rampFactor*t2*t24*t79*t102;
			f_udg[85*ng+i] += h*k_h*rampFactor*rz*t2*t24*t102;
			f_udg[86*ng+i] += h*k_h*rampFactor*ruz*t2*t24*t102;
			f_udg[87*ng+i] += h*k_h*rampFactor*rvz*t2*t24*t102;
			f_udg[88*ng+i] += h*k_h*rampFactor*rwz*t2*t24*t102;
			f_udg[89*ng+i] += h*k_h*rampFactor*t2*t24*t90*t102;
			f_udg[90*ng+i] += -h*k_h*rampFactor*rx*t2*t24*t106;
			f_udg[91*ng+i] += t103-h*k_h*rampFactor*rux*t2*t24*t106;
			f_udg[92*ng+i] += -h*k_h*rampFactor*rvx*t2*t24*t106;
			f_udg[93*ng+i] += -h*k_h*rampFactor*rwx*t2*t24*t106;
			f_udg[94*ng+i] += -h*k_h*rampFactor*t2*t24*t68*t106-h*k_h*rampFactor*ru*t2*t3*t24*t41*t51*2.0;
			f_udg[95*ng+i] += -h*k_h*rampFactor*ry*t2*t24*t106;
			f_udg[96*ng+i] += -h*k_h*rampFactor*ruy*t2*t24*t106;
			f_udg[97*ng+i] += -h*k_h*rampFactor*rvy*t2*t24*t106;
			f_udg[98*ng+i] += -h*k_h*rampFactor*rwy*t2*t24*t106;
			f_udg[99*ng+i] += -h*k_h*rampFactor*t2*t24*t79*t106;
			f_udg[100*ng+i] += -h*k_h*rampFactor*rz*t2*t24*t106;
			f_udg[101*ng+i] += -h*k_h*rampFactor*ruz*t2*t24*t106;
			f_udg[102*ng+i] += -h*k_h*rampFactor*rvz*t2*t24*t106;
			f_udg[103*ng+i] += -h*k_h*rampFactor*rwz*t2*t24*t106;
			f_udg[104*ng+i] += -h*k_h*rampFactor*t2*t24*t90*t106;
			f_udg[105*ng+i] += 0.0;
			f_udg[106*ng+i] += 0.0;
			f_udg[107*ng+i] += t103;
			f_udg[108*ng+i] += 0.0;
			f_udg[109*ng+i] += h*k_h*rampFactor*rv*t2*t3*t24*t41*t51*-2.0;
			f_udg[110*ng+i] += 0.0;
			f_udg[111*ng+i] += 0.0;
			f_udg[112*ng+i] += 0.0;
			f_udg[113*ng+i] += 0.0;
			f_udg[114*ng+i] += 0.0;
			f_udg[115*ng+i] += 0.0;
			f_udg[116*ng+i] += 0.0;
			f_udg[117*ng+i] += 0.0;
			f_udg[118*ng+i] += 0.0;
			f_udg[119*ng+i] += 0.0;
			f_udg[120*ng+i] += 0.0;
			f_udg[121*ng+i] += 0.0;
			f_udg[122*ng+i] += 0.0;
			f_udg[123*ng+i] += t103;
			f_udg[124*ng+i] += h*k_h*rampFactor*rw*t2*t3*t24*t41*t51*-2.0;
			f_udg[125*ng+i] += 0.0;
			f_udg[126*ng+i] += 0.0;
			f_udg[127*ng+i] += 0.0;
			f_udg[128*ng+i] += 0.0;
			f_udg[129*ng+i] += 0.0;
			f_udg[130*ng+i] += 0.0;
			f_udg[131*ng+i] += 0.0;
			f_udg[132*ng+i] += 0.0;
			f_udg[133*ng+i] += 0.0;
			f_udg[134*ng+i] += 0.0;
			f_udg[135*ng+i] += 0.0;
			f_udg[136*ng+i] += 0.0;
			f_udg[137*ng+i] += 0.0;
			f_udg[138*ng+i] += 0.0;
			f_udg[139*ng+i] += t115;
			f_udg[140*ng+i] += 0.0;
			f_udg[141*ng+i] += 0.0;
			f_udg[142*ng+i] += 0.0;
			f_udg[143*ng+i] += 0.0;
			f_udg[144*ng+i] += 0.0;
			f_udg[145*ng+i] += 0.0;
			f_udg[146*ng+i] += 0.0;
			f_udg[147*ng+i] += 0.0;
			f_udg[148*ng+i] += 0.0;
			f_udg[149*ng+i] += 0.0;
			f_udg[150*ng+i] += h*k_h*rampFactor*rx*t2*t24*t109;
			f_udg[151*ng+i] += h*k_h*rampFactor*rux*t2*t24*t109;
			f_udg[152*ng+i] += h*k_h*rampFactor*rvx*t2*t24*t109;
			f_udg[153*ng+i] += h*k_h*rampFactor*rwx*t2*t24*t109;
			f_udg[154*ng+i] += h*k_h*rampFactor*t2*t24*t68*t109;
			f_udg[155*ng+i] += t103+h*k_h*rampFactor*ry*t2*t24*t109;
			f_udg[156*ng+i] += h*k_h*rampFactor*ruy*t2*t24*t109;
			f_udg[157*ng+i] += h*k_h*rampFactor*rvy*t2*t24*t109;
			f_udg[158*ng+i] += h*k_h*rampFactor*rwy*t2*t24*t109;
			f_udg[159*ng+i] += t114+h*k_h*rampFactor*t2*t24*t79*t109;
			f_udg[160*ng+i] += h*k_h*rampFactor*rz*t2*t24*t109;
			f_udg[161*ng+i] += h*k_h*rampFactor*ruz*t2*t24*t109;
			f_udg[162*ng+i] += h*k_h*rampFactor*rvz*t2*t24*t109;
			f_udg[163*ng+i] += h*k_h*rampFactor*rwz*t2*t24*t109;
			f_udg[164*ng+i] += h*k_h*rampFactor*t2*t24*t90*t109;
			f_udg[165*ng+i] += 0.0;
			f_udg[166*ng+i] += 0.0;
			f_udg[167*ng+i] += 0.0;
			f_udg[168*ng+i] += 0.0;
			f_udg[169*ng+i] += 0.0;
			f_udg[170*ng+i] += 0.0;
			f_udg[171*ng+i] += t103;
			f_udg[172*ng+i] += 0.0;
			f_udg[173*ng+i] += 0.0;
			f_udg[174*ng+i] += h*k_h*rampFactor*ru*t2*t3*t24*t41*t51*-2.0;
			f_udg[175*ng+i] += 0.0;
			f_udg[176*ng+i] += 0.0;
			f_udg[177*ng+i] += 0.0;
			f_udg[178*ng+i] += 0.0;
			f_udg[179*ng+i] += 0.0;
			f_udg[180*ng+i] += -h*k_h*rampFactor*rx*t2*t24*t106;
			f_udg[181*ng+i] += -h*k_h*rampFactor*rux*t2*t24*t106;
			f_udg[182*ng+i] += -h*k_h*rampFactor*rvx*t2*t24*t106;
			f_udg[183*ng+i] += -h*k_h*rampFactor*rwx*t2*t24*t106;
			f_udg[184*ng+i] += -h*k_h*rampFactor*t2*t24*t68*t106;
			f_udg[185*ng+i] += -h*k_h*rampFactor*ry*t2*t24*t106;
			f_udg[186*ng+i] += -h*k_h*rampFactor*ruy*t2*t24*t106;
			f_udg[187*ng+i] += t103-h*k_h*rampFactor*rvy*t2*t24*t106;
			f_udg[188*ng+i] += -h*k_h*rampFactor*rwy*t2*t24*t106;
			f_udg[189*ng+i] += -h*k_h*rampFactor*t2*t24*t79*t106-h*k_h*rampFactor*rv*t2*t3*t24*t41*t51*2.0;
			f_udg[190*ng+i] += -h*k_h*rampFactor*rz*t2*t24*t106;
			f_udg[191*ng+i] += -h*k_h*rampFactor*ruz*t2*t24*t106;
			f_udg[192*ng+i] += -h*k_h*rampFactor*rvz*t2*t24*t106;
			f_udg[193*ng+i] += -h*k_h*rampFactor*rwz*t2*t24*t106;
			f_udg[194*ng+i] += -h*k_h*rampFactor*t2*t24*t90*t106;
			f_udg[195*ng+i] += 0.0;
			f_udg[196*ng+i] += 0.0;
			f_udg[197*ng+i] += 0.0;
			f_udg[198*ng+i] += 0.0;
			f_udg[199*ng+i] += 0.0;
			f_udg[200*ng+i] += 0.0;
			f_udg[201*ng+i] += 0.0;
			f_udg[202*ng+i] += 0.0;
			f_udg[203*ng+i] += t103;
			f_udg[204*ng+i] += h*k_h*rampFactor*rw*t2*t3*t24*t41*t51*-2.0;
			f_udg[205*ng+i] += 0.0;
			f_udg[206*ng+i] += 0.0;
			f_udg[207*ng+i] += 0.0;
			f_udg[208*ng+i] += 0.0;
			f_udg[209*ng+i] += 0.0;
			f_udg[210*ng+i] += 0.0;
			f_udg[211*ng+i] += 0.0;
			f_udg[212*ng+i] += 0.0;
			f_udg[213*ng+i] += 0.0;
			f_udg[214*ng+i] += 0.0;
			f_udg[215*ng+i] += 0.0;
			f_udg[216*ng+i] += 0.0;
			f_udg[217*ng+i] += 0.0;
			f_udg[218*ng+i] += 0.0;
			f_udg[219*ng+i] += t115;
			f_udg[220*ng+i] += 0.0;
			f_udg[221*ng+i] += 0.0;
			f_udg[222*ng+i] += 0.0;
			f_udg[223*ng+i] += 0.0;
			f_udg[224*ng+i] += 0.0;
			f_udg[225*ng+i] += h*k_h*rampFactor*rx*t2*t24*t118;
			f_udg[226*ng+i] += h*k_h*rampFactor*rux*t2*t24*t118;
			f_udg[227*ng+i] += h*k_h*rampFactor*rvx*t2*t24*t118;
			f_udg[228*ng+i] += h*k_h*rampFactor*rwx*t2*t24*t118;
			f_udg[229*ng+i] += h*k_h*rampFactor*t2*t24*t68*t118;
			f_udg[230*ng+i] += h*k_h*rampFactor*ry*t2*t24*t118;
			f_udg[231*ng+i] += h*k_h*rampFactor*ruy*t2*t24*t118;
			f_udg[232*ng+i] += h*k_h*rampFactor*rvy*t2*t24*t118;
			f_udg[233*ng+i] += h*k_h*rampFactor*rwy*t2*t24*t118;
			f_udg[234*ng+i] += h*k_h*rampFactor*t2*t24*t79*t118;
			f_udg[235*ng+i] += t103+h*k_h*rampFactor*rz*t2*t24*t118;
			f_udg[236*ng+i] += h*k_h*rampFactor*ruz*t2*t24*t118;
			f_udg[237*ng+i] += h*k_h*rampFactor*rvz*t2*t24*t118;
			f_udg[238*ng+i] += h*k_h*rampFactor*rwz*t2*t24*t118;
			f_udg[239*ng+i] += t114+h*k_h*rampFactor*t2*t24*t90*t118;
			f_udg[240*ng+i] += 0.0;
			f_udg[241*ng+i] += 0.0;
			f_udg[242*ng+i] += 0.0;
			f_udg[243*ng+i] += 0.0;
			f_udg[244*ng+i] += 0.0;
			f_udg[245*ng+i] += 0.0;
			f_udg[246*ng+i] += 0.0;
			f_udg[247*ng+i] += 0.0;
			f_udg[248*ng+i] += 0.0;
			f_udg[249*ng+i] += 0.0;
			f_udg[250*ng+i] += 0.0;
			f_udg[251*ng+i] += t103;
			f_udg[252*ng+i] += 0.0;
			f_udg[253*ng+i] += 0.0;
			f_udg[254*ng+i] += h*k_h*rampFactor*ru*t2*t3*t24*t41*t51*-2.0;
			f_udg[255*ng+i] += 0.0;
			f_udg[256*ng+i] += 0.0;
			f_udg[257*ng+i] += 0.0;
			f_udg[258*ng+i] += 0.0;
			f_udg[259*ng+i] += 0.0;
			f_udg[260*ng+i] += 0.0;
			f_udg[261*ng+i] += 0.0;
			f_udg[262*ng+i] += 0.0;
			f_udg[263*ng+i] += 0.0;
			f_udg[264*ng+i] += 0.0;
			f_udg[265*ng+i] += 0.0;
			f_udg[266*ng+i] += 0.0;
			f_udg[267*ng+i] += t103;
			f_udg[268*ng+i] += 0.0;
			f_udg[269*ng+i] += h*k_h*rampFactor*rv*t2*t3*t24*t41*t51*-2.0;
			f_udg[270*ng+i] += -h*k_h*rampFactor*rx*t2*t24*t106;
			f_udg[271*ng+i] += -h*k_h*rampFactor*rux*t2*t24*t106;
			f_udg[272*ng+i] += -h*k_h*rampFactor*rvx*t2*t24*t106;
			f_udg[273*ng+i] += -h*k_h*rampFactor*rwx*t2*t24*t106;
			f_udg[274*ng+i] += -h*k_h*rampFactor*t2*t24*t68*t106;
			f_udg[275*ng+i] += -h*k_h*rampFactor*ry*t2*t24*t106;
			f_udg[276*ng+i] += -h*k_h*rampFactor*ruy*t2*t24*t106;
			f_udg[277*ng+i] += -h*k_h*rampFactor*rvy*t2*t24*t106;
			f_udg[278*ng+i] += -h*k_h*rampFactor*rwy*t2*t24*t106;
			f_udg[279*ng+i] += -h*k_h*rampFactor*t2*t24*t79*t106;
			f_udg[280*ng+i] += -h*k_h*rampFactor*rz*t2*t24*t106;
			f_udg[281*ng+i] += -h*k_h*rampFactor*ruz*t2*t24*t106;
			f_udg[282*ng+i] += -h*k_h*rampFactor*rvz*t2*t24*t106;
			f_udg[283*ng+i] += t103-h*k_h*rampFactor*rwz*t2*t24*t106;
			f_udg[284*ng+i] += -h*k_h*rampFactor*t2*t24*t90*t106-h*k_h*rampFactor*rw*t2*t3*t24*t41*t51*2.0;
			f_udg[285*ng+i] += 0.0;
			f_udg[286*ng+i] += 0.0;
			f_udg[287*ng+i] += 0.0;
			f_udg[288*ng+i] += 0.0;
			f_udg[289*ng+i] += 0.0;
			f_udg[290*ng+i] += 0.0;
			f_udg[291*ng+i] += 0.0;
			f_udg[292*ng+i] += 0.0;
			f_udg[293*ng+i] += 0.0;
			f_udg[294*ng+i] += 0.0;
			f_udg[295*ng+i] += 0.0;
			f_udg[296*ng+i] += 0.0;
			f_udg[297*ng+i] += 0.0;
			f_udg[298*ng+i] += 0.0;
			f_udg[299*ng+i] += t115;
		}

	}
}

