void fbou_leuq2d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, 
        double *nl, double *ui, double *param, double time, int ib, int ng, int nc, int ncu, int nd, int ncd)
{
	double mu = param[0];
	double lambda = param[1];
	double tau = param[2];
    double uinf1 = ui[0];
    double uinf2 = ui[1];
    double t2, t3, t4, t5;
    
	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double q11 = udg[2*ng+i];
		double q21 = udg[3*ng+i];
		double q12 = udg[4*ng+i];
		double q22 = udg[5*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
        
        if (ib==1) {
            fh[0*ng+i] = -uh1+uinf1+x1;
            fh[1*ng+i] = -uh2+uinf2+x2;
            fh_udg[0*ng+i] = 0.0;
            fh_udg[1*ng+i] = 0.0;
            fh_udg[2*ng+i] = 0.0;
            fh_udg[3*ng+i] = 0.0;
            fh_udg[4*ng+i] = 0.0;
            fh_udg[5*ng+i] = 0.0;
            fh_udg[6*ng+i] = 0.0;
            fh_udg[7*ng+i] = 0.0;
            fh_udg[8*ng+i] = 0.0;
            fh_udg[9*ng+i] = 0.0;
            fh_udg[10*ng+i] = 0.0;
            fh_udg[11*ng+i] = 0.0; 
            fh_uh[0*ng+i] = -1.0;
            fh_uh[3*ng+i] = -1.0;
        }
        else if (ib ==2) {
            t2 = q11+q22+2.0;
            t3 = lambda*t2;
            t4 = q12+q21;
            fh[0*ng+i] = -uinf1+nl1*(t3+mu*(q11+1.0)*2.0)+tau*(u1-uh1)+mu*nl2*t4;
            fh[1*ng+i] = -uinf2+nl2*(t3+mu*(q22+1.0)*2.0)+tau*(u2-uh2)+mu*nl1*t4;

            t2 = mu*nl2;
            t3 = mu*nl1;
            t4 = mu*2.0;
            t5 = lambda+t4;
            fh_udg[0*ng+i] = tau;
            fh_udg[1*ng+i] = 0.0;
            fh_udg[2*ng+i] = nl1*t5;
            fh_udg[3*ng+i] = t2;
            fh_udg[4*ng+i] = t2;
            fh_udg[5*ng+i] = lambda*nl1;
            fh_udg[6*ng+i] = 0.0;
            fh_udg[7*ng+i] = tau;
            fh_udg[8*ng+i] = lambda*nl2;
            fh_udg[9*ng+i] = t3;
            fh_udg[10*ng+i] = t3;
            fh_udg[11*ng+i] = nl2*t5;
            fh_uh[0*ng+i] = -tau;
            fh_uh[3*ng+i] = -tau;
        }
        else if (ib ==3) {
            fh[0*ng+i] = -uh1+uinf1+x1;
            fh[1*ng+i] = -uinf2+tau*(u2-uh2)+nl2*(lambda*(q11+q22+2.0)+mu*(q22+1.0)*2.0)+mu*nl1*(q12+q21);

            t2 = mu*nl1;
            fh_udg[0*ng+i] = 0.0;
            fh_udg[1*ng+i] = 0.0;
            fh_udg[2*ng+i] = 0.0;
            fh_udg[3*ng+i] = 0.0;
            fh_udg[4*ng+i] = 0.0;
            fh_udg[5*ng+i] = 0.0;
            fh_udg[6*ng+i] = 0.0;
            fh_udg[7*ng+i] = tau;
            fh_udg[8*ng+i] = lambda*nl2;
            fh_udg[9*ng+i] = t2;
            fh_udg[10*ng+i] = t2;
            fh_udg[11*ng+i] = nl2*(lambda+mu*2.0);

            fh_uh[0*ng+i] = -1.0
            fh_uh[3*ng+i] = -tau
        }
        else if (ib ==4) {
            fh[0*ng+i] = -uinf1+tau*(u1-uh1)+nl1*(lambda*(q11+q22+2.0)+mu*(q11+1.0)*2.0)+mu*nl2*(q12+q21);
            fh[1*ng+i] = -uh2+uinf2+x2;

            t2 = mu*nl2;
            fh_udg[0*ng+i] = tau;
            fh_udg[1*ng+i] = 0.0;
            fh_udg[2*ng+i] = nl1*(lambda+mu*2.0);
            fh_udg[3*ng+i] = t2;
            fh_udg[4*ng+i] = t2;
            fh_udg[5*ng+i] = lambda*nl1;
            fh_udg[6*ng+i] = 0.0;
            fh_udg[7*ng+i] = 0.0;
            fh_udg[8*ng+i] = 0.0;
            fh_udg[9*ng+i] = 0.0;
            fh_udg[10*ng+i] = 0.0;
            fh_udg[11*ng+i] = 0.0;

            fh_uh[0*ng+i] = -tau;
            fh_uh[3*ng+i] = -1.0;
        }
	}
}

