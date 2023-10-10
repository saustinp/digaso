#ifndef __ELEMFACEGEOM
#define __ELEMFACEGEOM

void elementgeom(double *p, double *J, double *Xx, double *jac, double *pn, double *shapvt, Int *ndims, Int numPoints)
{
//     numPoints: ngv = ndims[13];
    
    Int i, iA, iB, inc = 1;
    char chn = 'N';
    double one = 1.0, zero = 0.0;
    
    /* Get dimensions */
    Int nd, npv, ncd;
    nd = ndims[0];
    ncd = ndims[1];
    npv = ndims[9];
    
    DGEMM(&chn, &chn, &numPoints, &ncd, &npv, &one, shapvt, &numPoints, pn,
          &npv, &zero, p, &numPoints);
    /* shapvt: numPoints / npv / (1+nd) (z.r) (we only use numPoints / npv here) */
    /* pn: npv / ncd */
    /* p: numPoints / ncd */
    
    /* Compute Jacobian matrix at desired points */
    for (i=0; i<nd; i++) {
        iA = (i+1)*numPoints*npv;
        iB = i*numPoints*nd;
        DGEMM(&chn, &chn, &numPoints, &nd, &npv, &one, &shapvt[iA], &numPoints, pn,
              &npv, &zero, &J[iB], &numPoints);
    }
    /* shapvt: numPoints / npv / (1+nd) (z.r) (we only use numPoints / npv / nd (z.r) here) */
    /* pn: npv / ncd (we only use npv / nd (z.x) here) */
    /* J: numPoints / nd (z.x) / nd (z.r) */

    double *J11, *J12, *J13, *J21, *J22, *J23, *J31, *J32, *J33;
    double *Xx11, *Xx12, *Xx13, *Xx21, *Xx22, *Xx23, *Xx31, *Xx32, *Xx33;

    /* Compute determinant and inverse times Jacobian of mapping from master element to undeformed element evaluated at desired points */
    if (nd==1) {
        for (i=0; i<numPoints; i++) {
            jac[i] = J[i];
            Xx[i] = 1.0;
        }
    }
    else if (nd==2) {
        J11 = &J[0*numPoints];
        J12 = &J[1*numPoints];
        J21 = &J[2*numPoints];
        J22 = &J[3*numPoints];
        Xx11 = &Xx[0*numPoints];
        Xx21 = &Xx[1*numPoints];
        Xx12 = &Xx[2*numPoints];
        Xx22 = &Xx[3*numPoints];
        for (i=0; i<numPoints; i++) {
            jac[i] = J11[i]*J22[i] - J12[i]*J21[i];
            Xx11[i] = J22[i];
            Xx21[i] = -J21[i];
            Xx12[i] = -J12[i];
            Xx22[i] = J11[i];
        }
    }
    else if (nd==3) {
        J11 = &J[0*numPoints];
        J12 = &J[1*numPoints];
        J13 = &J[2*numPoints];
        J21 = &J[3*numPoints];
        J22 = &J[4*numPoints];
        J23 = &J[5*numPoints];
        J31 = &J[6*numPoints];
        J32 = &J[7*numPoints];
        J33 = &J[8*numPoints];
        Xx11 = &Xx[0*numPoints];
        Xx21 = &Xx[1*numPoints];
        Xx31 = &Xx[2*numPoints];
        Xx12 = &Xx[3*numPoints];
        Xx22 = &Xx[4*numPoints];
        Xx32 = &Xx[5*numPoints];
        Xx13 = &Xx[6*numPoints];
        Xx23 = &Xx[7*numPoints];
        Xx33 = &Xx[8*numPoints];
        for (i=0; i<numPoints; i++) {
            jac[i] = J11[i]*J22[i]*J33[i] - J11[i]*J32[i]*J23[i] +
                     J21[i]*J32[i]*J13[i] - J21[i]*J12[i]*J33[i] +
                     J31[i]*J12[i]*J23[i] - J31[i]*J22[i]*J13[i];
            Xx11[i] = J22[i]*J33[i] - J23[i]*J32[i];
            Xx21[i] = J23[i]*J31[i] - J21[i]*J33[i];
            Xx31[i] = J21[i]*J32[i] - J22[i]*J31[i];
            Xx12[i] = J13[i]*J32[i] - J12[i]*J33[i];
            Xx22[i] = J11[i]*J33[i] - J13[i]*J31[i];
            Xx32[i] = J12[i]*J31[i] - J11[i]*J32[i];
            Xx13[i] = J12[i]*J23[i] - J13[i]*J22[i];
            Xx23[i] = J13[i]*J21[i] - J11[i]*J23[i];
            Xx33[i] = J11[i]*J22[i] - J12[i]*J21[i];
        }
    }
    /* jac: numPoints */
    /* Xx: numPoints / nd (z.r) / nd (z.x) */
}

void facegeom(double *p, double *J, double *Xx, double* jac, double *nl, 
        double *pf, double *shapft, Int *ndims, Int numPoints)
{    
    // Hint: pf needs to be precomputed before
    Int inc = 1, i, j, k, iA, iB, is, na, nb;
    char chn = 'N';
    double one = 1.0, zero = 0.0;

    /* Get dimensions */
    Int nd, npv, npf, nfe, ncd, ndf;
    nd = ndims[0];
    ncd = ndims[1];
    nfe = ndims[2];
    npv = ndims[9];
    npf = ndims[10];
    ndf = npf*nfe;

    /* Compute nodal fields at desired points */
    na = nfe*ncd;
    DGEMM(&chn, &chn, &numPoints, &na, &npf, &one, shapft, &numPoints, pf,
          &npf, &zero, p, &numPoints);
    /* shapft: numPoints / npf / nd (z.r) (only numPoints / npf is used here) */
    /* pf: npf / nfe / ncd */
    /* p: numPoints / nfe / ncd */

    /* Compute Jacobian matrix at desired points */
    nb = nfe*nd;
    for (i=0; i<nd-1; i++) {
        iA = (i+1)*numPoints*npf;
        iB = i*numPoints*nfe*nd;
        DGEMM(&chn, &chn, &numPoints, &nb, &npf, &one, &shapft[iA], &numPoints, pf,
              &npf, &zero, &J[iB], &numPoints);
    }
    /* shapft: numPoints / npf / nd (z.r) (only numPoints / npf / nd-1 is used here) */
    /* pf: npf / nfe / ncd (only npf / nfe / nd (z.x) is used here) */
    /* J: numPoints / nfe / nd (z.x) / nd-1 (z.r) */
        
    double *J11, *J12, *J21, *J22, *J31, *J32;

    /* Compute determinant and normal vector at desired points */
    na = numPoints*nfe;
    if (nd==1) {
        jac[0] = 1.0;
        nl[0] = 1.0;
        nl[1] = -1.0;
    }
    else if (nd==2) {
        for (i=0; i<na; i++) {
            j = i+na;
            jac[i] = sqrt(J[i]*J[i] + J[j]*J[j]);
            nl[i] = J[j]/jac[i];
            nl[j] = -J[i]/jac[i];
        }
    }
    else if (nd==3) {
        J11 = &J[0*na];
        J21 = &J[1*na];
        J31 = &J[2*na];
        J12 = &J[3*na];
        J22 = &J[4*na];
        J32 = &J[5*na];
        for (i=0; i<na; i++) {
            j = i+na;
            k = i+2*na;
            nl[i] = J21[i]*J32[i] - J31[i]*J22[i];
            nl[j] = J31[i]*J12[i] - J11[i]*J32[i];
            nl[k] = J11[i]*J22[i] - J21[i]*J12[i];
            jac[i] = sqrt(nl[i]*nl[i] + nl[j]*nl[j] + nl[k]*nl[k]);
            nl[i] = nl[i]/jac[i];
            nl[j] = nl[j]/jac[i];
            nl[k] = nl[k]/jac[i];
        }
    }
    /* J: numPoints / nfe / nd (z.x) / nd-1 (z.r). Except in 1D, that is length 2. */
    /* Nfg: numPoints / nfe / nd (z.x) */
    /* jac: numPoints / nfe. Except in 1D, that is length 1. */
    /* p: numPoints / nfe / ncd */
}

#endif
