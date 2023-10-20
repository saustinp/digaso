void sourceonly_axispoisson2d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	// Physics parameters
	double l_ref = param[1];
	double E_ref = param[3];
	double e_eps0 = param[4];
	double N0 = param[6];
	double z0 = param[7];
	double sigma0 = param[8];

	double N0_tilde = N0*(l_ref*l_ref*l_ref);
	double z0_tilde = z0/l_ref;
	double sigma0_tilde = sigma0/l_ref;

	int i;
	for (i = 0; i <ng*ncu; i++) {
		double r_tilde = pg[0*ng + i];
		double z_tilde = pg[1*ng + i];
		s[i] = r_tilde*(e_eps0/(E_ref*l_ref*l_ref))*N0_tilde*exp(-((z_tilde - z0_tilde)*(z_tilde - z0_tilde) + r_tilde*r_tilde)/(sigma0_tilde*sigma0_tilde));
		
	}
}

