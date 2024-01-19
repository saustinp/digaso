
void fbou2d(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nc, int nch, int nd, int ncd)
{            
    
    int    i;                   
    double tau = param[9];   
	double r1=1.0;       // Mod for axisymmetric

    if (ib==1) { // Electrode: species set to homogeneous nuemann and potential set to dirichlet
		
		// Zero all arrays
		for (i=0; i<ng*nch; i++)
			fh[i] = 0.0;
		for (i=0; i<ng*nch*nc; i++)
			fh_u[i] = 0.0;
		for (i=0; i<ng*nch*nch; i++)
			fh_uh[i] = 0.0;

		fhat2d(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd);

// 		for (i=0; i<ng; i++)
// 			fh[2*ng + i] = r1*tau*(ui[i]-uhg[2*ng+i]);
    
// 		for (int k=0; k<nc; k++){
// 			for (int j=0; j<nch; j++){
// 				for (int i=0; i<ng; i++) {
// 					if (j == 2)
// 						fh_u[ng*nch*k + ng*j + i] = 0.0;    // Potential dirichlet
// 				}
// 			}
// 		}
// 
// 		for (int k=0; k<nch; k++){
// 			for (int j=0; j<nch; j++){
// 				for (int i=0; i<ng; i++) {
// 					if (j == 2)
// 						fh_uh[ng*nch*k + ng*j + i] = -r1*tau;    // Potential dirichlet
// 				}
// 			}
//     }

		for (i=0; i<ng; i++)
			fh[2*ng + i] = r1*tau*(ui[2]-uhg[2*ng+i]);

		for (int k=0; k<nc; k++){
			for (int j=0; j<nch; j++){
				for (int i=0; i<ng; i++) {
					if (j == (nch-1))
						fh_u[ng*nch*k + ng*j + i] = 0.0;    // Potential dirichlet
				}
			}
		}

		for (int k=0; k<nch; k++){
			for (int j=0; j<nch; j++){
				for (int i=0; i<ng; i++) {
					if (j == (nch-1)) {
            if (j==k)
              fh_uh[ng*nch*k + ng*j + i] = -r1*tau;    // Potential dirichlet
            else
              fh_uh[ng*nch*k + ng*j + i] = 0.0;
          }
				}
			}
    }    
	}

    else if (ib==2) { // Right farfield -- species + potential all have homogeneous neumann
        fhat2d(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd);    

    }  

    else if (ib==3) { // Symmetry boundary: extrapolate m=u or u=uhat
        // Zero all arrays
        for (i=0; i<ng*nch; i++)
            fh[i] = 0.0;
        for (i=0; i<ng*nch*nc; i++)
            fh_u[i] = 0.0;
        for (i=0; i<ng*nch*nch; i++)
            fh_uh[i] = 0.0;

        for (i=0; i<ng; i++) {

            // fh(:,1), fh(:,2), fh(:,3)
            fh[0*ng+i] = r1*tau*(udg[0*ng+i]-uhg[0*ng+i]);        // Check that uhg is correct here?
            fh[1*ng+i] = r1*tau*(udg[1*ng+i]-uhg[1*ng+i]);
            fh[2*ng+i] = r1*tau*(udg[2*ng+i]-uhg[2*ng+i]);
            
            // fh_udg(:,1,1), fh_udg(:,2,2), fh_udg(:,3,3)
            fh_u[0*ng*nch+0*ng+i] = r1*tau;
            fh_u[1*ng*nch+1*ng+i] = r1*tau;
            fh_u[2*ng*nch+2*ng+i] = r1*tau;

            // fh_uh(:,1,1), fh_uh(:,2,2), fh_uh(:,3,3)
            fh_uh[0*ng*nch+0*ng+i] = -r1*tau;
            fh_uh[1*ng*nch+1*ng+i] = -r1*tau;
            fh_uh[2*ng*nch+2*ng+i] = -r1*tau;
        }
	}

}