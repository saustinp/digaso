#ifndef __STABILIZATIONTENSOR
#define __STABILIZATIONTENSOR

#include "../../utilities/UH2uh.cpp"
#include "../../utilities/NL2nl.cpp"

// Written by: C. Nguyen & P. Fernandez

void getConstantStabilizationTensor(double* tau, double* tau_uh, double* param, int numPoints, int nch, int computeJacobian)
{
    int i, j, k;
    int sz, sz2 = numPoints * nch;

    double tauValue = param[5];

    for (k=0; k<nch; k++)
        for (j=0; j<nch; j++)
            for (i = 0; i < numPoints; i++)
                tau[i + j * numPoints + k * sz2] = (j==k) ? tauValue : 0.0;
    // tau: numPoints / nch / nch

    if (computeJacobian == 1) {
        sz = numPoints * nch * nch * nch;
        for (i=0; i<sz; i++)
            tau_uh[i] = 0.0;
    }
    // tau_uh: numPoints / nch / nch / nch
}


void getLaxFriedrichStabilizationTensor2d(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    // This functions works for compressible Euler, NS and RANS (any turbulence model).
    // If tauValue <= 0, the maximum eigenvalue is used to stabilize the turbulence model equation(s). Otherwise, tauValue is used.

    double gam, gam1, epslm, tauValue, Dc2Reg_Dc2;
    double nx, ny;
    double r, ru, rv, rE;
    double r1, r1m1, r1m2, r1m3, r1m4;
    double uv, uvm1, uvm2, uvm3, uvm4;
    double vv, vvm1, vvm2, vvm3, vvm4;
    double Vgx, Vgy, Vgn;
    double af, afm1, afm2, afm3, afm4;
    double p, pm1, pm2, pm3, pm4;
    double c2, c2m1, c2m2, c2m3, c2m4;
    double c, cm1, cm2, cm3, cm4;
    double un, unm1, unm2, unm3, unm4;
    double lamMax, lamMaxm1, lamMaxm2, lamMaxm3, lamMaxm4;
    double DlamMax1Dx;

    gam   = param[0];
    gam1 = gam - 1.0;
    epslm = 0.2;//param[1];
    tauValue = param[5];
    double pi = 3.141592653589793;
    double b = 100.0;
    double minTau = 0.01;
    
    int i, j, k, m;
    int sz2 = nch * numPoints;

    double stabTurbEq;
    double DlamMax[nch];
    double DstabTurbEq[nch];

    Int ALEflag = app.ALEflag;

    for (i=0; i<numPoints; i++) {
        r   = uh[0*numPoints+i]; // 1.0000;
        ru  = uh[1*numPoints+i]; // 0.9994;
        rv  = uh[2*numPoints+i]; // 0.0349;
        rE  = uh[3*numPoints+i]; // 45.1429;
        
        nx = nl[0 * numPoints + i];
        ny = nl[1 * numPoints + i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
        }

        r1 = 1.0 / r;
        uv = ru * r1;
        vv = rv * r1;
        af = 0.5 * (uv * uv + vv * vv);
        p = gam1 * (rE - r * af);
        c2 = gam * p * r1;
        Dc2Reg_Dc2 = atan(b*c2)/pi + 1.0/2.0 + (b*c2) / ((1+b*c2*b*c2) * pi);
        c2 = max ( c2*(atan(b*c2)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 , 1.0e-6);
        c = sqrt(c2);
        un = uv * nx + vv * ny;

        if (computeJacobian == 1) {
            r1m1 = -1.0 / (r * r);
            r1m2 = 0.0;
            r1m3 = 0.0;
            r1m4 = 0.0;

            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;
            uvm3 = ru * r1m3;
            uvm4 = ru * r1m4;

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;
            vvm3 = r1 + rv * r1m3;
            vvm4 = rv * r1m4;

            afm1 = uv * uvm1 + vv * vvm1;
            afm2 = uv * uvm2 + vv * vvm2;
            afm3 = uv * uvm3 + vv * vvm3;
            afm4 = uv * uvm4 + vv * vvm4;

            pm1 = gam1 * (-af - r * afm1);
            pm2 = gam1 * (-r * afm2);
            pm3 = gam1 * (-r * afm3);
            pm4 = gam1 * (1.0 - r * afm4);

            c2m1 = gam * (pm1 * r1 + p * r1m1);
            c2m2 = gam * (pm2 * r1 + p * r1m2);
            c2m3 = gam * (pm3 * r1 + p * r1m3);
            c2m4 = gam * (pm4 * r1 + p * r1m4);
            
            c2m1 *= Dc2Reg_Dc2;
            c2m2 *= Dc2Reg_Dc2;
            c2m3 *= Dc2Reg_Dc2;
            c2m4 *= Dc2Reg_Dc2;

            cm1 = 0.5 * c2m1 / c;
            cm2 = 0.5 * c2m2 / c;
            cm3 = 0.5 * c2m3 / c;
            cm4 = 0.5 * c2m4 / c;

            unm1 = uvm1 * nx + vvm1 * ny;
            unm2 = uvm2 * nx + vvm2 * ny;
            unm3 = uvm3 * nx + vvm3 * ny;
            unm4 = uvm4 * nx + vvm4 * ny;
        }
        
        Vgn = Vgx * nx + Vgy * ny;
        
        if (un - Vgn >= 0) {
            lamMax = (un-Vgn) + c - minTau;
            DlamMax1Dx = 2*(atan(b*lamMax)/pi + 1.0/2.0 + (b*lamMax) / ((1+b*lamMax*b*lamMax) * pi)) - 1.0;
            lamMax = 2*(lamMax*(atan(b*lamMax)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0) - lamMax + minTau;
            if (computeJacobian == 1) {
                lamMaxm1 = DlamMax1Dx * (unm1 + cm1);
                lamMaxm2 = DlamMax1Dx * (unm2 + cm2);
                lamMaxm3 = DlamMax1Dx * (unm3 + cm3);
                lamMaxm4 = DlamMax1Dx * (unm4 + cm4);
            }
        }
        else {
            lamMax = - (un-Vgn) + c - minTau;
            DlamMax1Dx = 2*(atan(b*lamMax)/pi + 1.0/2.0 + (b*lamMax) / ((1+b*lamMax*b*lamMax) * pi)) - 1.0;
            lamMax = 2*(lamMax*(atan(b*lamMax)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0) - lamMax + minTau;
            if (computeJacobian == 1) {
                lamMaxm1 =  DlamMax1Dx * (-unm1 + cm1);
                lamMaxm2 =  DlamMax1Dx * (-unm2 + cm2);
                lamMaxm3 =  DlamMax1Dx * (-unm3 + cm3);
                lamMaxm4 =  DlamMax1Dx * (-unm4 + cm4);
            }
        }

       if (epslm > 0) {
           lamMax = 0.5 * (lamMax * lamMax / (epslm * c) + epslm * c);
           lamMaxm1 = lamMax * lamMaxm1 / (epslm*c) + 0.5 * (1 - lamMax * lamMax / (epslm * c * epslm * c)) * (epslm * cm1);
           lamMaxm2 = lamMax * lamMaxm2 / (epslm*c) + 0.5 * (1 - lamMax * lamMax / (epslm * c * epslm * c)) * (epslm * cm2);
           lamMaxm3 = lamMax * lamMaxm3 / (epslm*c) + 0.5 * (1 - lamMax * lamMax / (epslm * c * epslm * c)) * (epslm * cm3);
           lamMaxm4 = lamMax * lamMaxm4 / (epslm*c) + 0.5 * (1 - lamMax * lamMax / (epslm * c * epslm * c)) * (epslm * cm4);
       }

        for (k = 0; k < 4; k++)
            for (j = 0; j < nch; j++)
                tau[i + j * numPoints + k * sz2] = (j == k) ? lamMax : 0.0;
        
        if (tauValue > 0) {
            stabTurbEq = tauValue;
            if (computeJacobian == 1) {
                for (j = 0; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }
        else {
            stabTurbEq = lamMax;
            if (computeJacobian == 1) {
                DstabTurbEq[0] = lamMaxm1;
                DstabTurbEq[1] = lamMaxm2;
                DstabTurbEq[2] = lamMaxm3;
                DstabTurbEq[3] = lamMaxm4;
                for (j = 4; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }

        for (k = 4; k < nch; k++)
            for (j = 0; j < nch; j++)
                tau[i + j * numPoints + k * sz2] = (j == k) ? stabTurbEq : 0.0;
        // tau: numPoints / nch / nch

        if (computeJacobian == 1) {
            int sz3 = nch * sz2;
            DlamMax[0] = lamMaxm1;
            DlamMax[1] = lamMaxm2;
            DlamMax[2] = lamMaxm3;
            DlamMax[3] = lamMaxm4;
            for (j=4; j<nch; j++)
                DlamMax[j] = 0.0;

            for (m = 0; m < nch; m++)
                for (k = 0; k < 4; k++)
                    for (j = 0; j < nch; j++)
                        tau_uh[i + j * numPoints + k * sz2 + m * sz3] = (j == k) ? DlamMax[m] : 0.0;

            for (m = 0; m < nch; m++)
                for (k = 4; k < nch; k++)
                    for (j = 0; j < nch; j++)
                        tau_uh[i + j * numPoints + k * sz2 + m * sz3] = (j == k) ? DstabTurbEq[m] : 0.0;
        }
        // tau_uh: numPoints / nch / nch / nch
    }
    
//     for (i=0; i<numPoints*nch*nch*nch; i++)
//         tau_uh[i] = 0.0;
}


void getLaxFriedrichStabilizationTensor3d(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    // This functions works for compressible Euler, NS and RANS (any turbulence model).
    // If tauValue <= 0, the maximum eigenvalue is used to stabilize the turbulence model equation(s). Otherwise, tauValue is used.

    double gam, gam1, epslm, tauValue, Dc2Reg_Dc2;
    double nx, ny, nz;
    double r, ru, rv, rw, rE;
    double r1, r1m1, r1m2, r1m3, r1m4, r1m5;
    double uv, uvm1, uvm2, uvm3, uvm4, uvm5;
    double vv, vvm1, vvm2, vvm3, vvm4, vvm5;
    double wv, wvm1, wvm2, wvm3, wvm4, wvm5;
    double Vgx, Vgy, Vgz, Vgn;
    double af, afm1, afm2, afm3, afm4, afm5;
    double p, pm1, pm2, pm3, pm4, pm5;
    double c2, c2m1, c2m2, c2m3, c2m4, c2m5;
    double c, cm1, cm2, cm3, cm4, cm5;
    double un, unm1, unm2, unm3, unm4, unm5;
    double lamMax, lamMaxm1, lamMaxm2, lamMaxm3, lamMaxm4, lamMaxm5;
    double DlamMax1Dx;

    gam = param[0];
    gam1 = gam - 1.0;
    epslm = param[1];
    tauValue = param[5];
    double pi = 3.141592653589793;
    double b = 100.0;
    double minTau = 0.1;

    int i, j, k, m;
    int sz2 = nch * numPoints;

    double stabTurbEq;
    double DlamMax[nch];
    double DstabTurbEq[nch];

    Int ALEflag = app.ALEflag;

    for (i = 0; i < numPoints; i++) {
        nx = nl[0 * numPoints + i];
        ny = nl[1 * numPoints + i];
        nz = nl[2 * numPoints + i];

        r = uh[0 * numPoints + i];
        ru = uh[1 * numPoints + i];
        rv = uh[2 * numPoints + i];
        rw = uh[3 * numPoints + i];
        rE = uh[4 * numPoints + i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
            Vgz = pg[(2 * nd + 2) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
            Vgz = 0.0;
        }

        r1 = 1.0 / r;
        uv = ru * r1;
        vv = rv * r1;
        wv = rw * r1;
        af = 0.5 * (uv * uv + vv * vv + wv * wv);
        p = gam1 * (rE - r * af);
        c2 = gam * p * r1;
        Dc2Reg_Dc2 = atan(b*c2)/pi + 1.0/2.0 + (b*c2) / ((1+b*c2*b*c2) * pi);
        c2 = max ( c2*(atan(b*c2)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 , 1.0e-6);
        c = sqrt(c2);
        un = uv * nx + vv * ny + wv * nz;

        if (computeJacobian == 1) {
            r1m1 = -1.0 / (r * r);
            r1m2 = 0.0;
            r1m3 = 0.0;
            r1m4 = 0.0;
            r1m5 = 0.0;

            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;
            uvm3 = ru * r1m3;
            uvm4 = ru * r1m4;
            uvm5 = ru * r1m5;

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;
            vvm3 = r1 + rv * r1m3;
            vvm4 = rv * r1m4;
            vvm5 = rv * r1m5;

            wvm1 = rw * r1m1;
            wvm2 = rw * r1m2;
            wvm3 = rw * r1m3;
            wvm4 = r1 + rw * r1m4;
            wvm5 = rw * r1m5;

            afm1 = uv * uvm1 + vv * vvm1 + wv * wvm1;
            afm2 = uv * uvm2 + vv * vvm2 + wv * wvm2;
            afm3 = uv * uvm3 + vv * vvm3 + wv * wvm3;
            afm4 = uv * uvm4 + vv * vvm4 + wv * wvm4;
            afm5 = uv * uvm5 + vv * vvm5 + wv * wvm5;

            pm1 = gam1 * (-af - r * afm1);
            pm2 = gam1 * (-r * afm2);
            pm3 = gam1 * (-r * afm3);
            pm4 = gam1 * (-r * afm4);
            pm5 = gam1 * (1.0 - r * afm5);

            c2m1 = gam * (pm1 * r1 + p * r1m1);
            c2m2 = gam * (pm2 * r1 + p * r1m2);
            c2m3 = gam * (pm3 * r1 + p * r1m3);
            c2m4 = gam * (pm4 * r1 + p * r1m4);
            c2m5 = gam * (pm5 * r1 + p * r1m5);
            
            c2m1 *= Dc2Reg_Dc2;
            c2m2 *= Dc2Reg_Dc2;
            c2m3 *= Dc2Reg_Dc2;
            c2m4 *= Dc2Reg_Dc2;
            c2m5 *= Dc2Reg_Dc2;

            cm1 = 0.5 * c2m1 / c;
            cm2 = 0.5 * c2m2 / c;
            cm3 = 0.5 * c2m3 / c;
            cm4 = 0.5 * c2m4 / c;
            cm5 = 0.5 * c2m5 / c;

            unm1 = uvm1 * nx + vvm1 * ny + wvm1 * nz;
            unm2 = uvm2 * nx + vvm2 * ny + wvm2 * nz;
            unm3 = uvm3 * nx + vvm3 * ny + wvm3 * nz;
            unm4 = uvm4 * nx + vvm4 * ny + wvm4 * nz;
            unm5 = uvm5 * nx + vvm5 * ny + wvm5 * nz;
        }
        
        Vgn = Vgx * nx + Vgy * ny + Vgz * nz;
        
        if (un - Vgn >= 0) {
            lamMax = (un-Vgn) + c - minTau;
            DlamMax1Dx = 2*(atan(b*lamMax)/pi + 1.0/2.0 + (b*lamMax) / ((1+b*lamMax*b*lamMax) * pi)) - 1.0;
            lamMax = 2*(lamMax*(atan(b*lamMax)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0) - lamMax + minTau;
            if (computeJacobian == 1) {
                lamMaxm1 = DlamMax1Dx * (unm1 + cm1);
                lamMaxm2 = DlamMax1Dx * (unm2 + cm2);
                lamMaxm3 = DlamMax1Dx * (unm3 + cm3);
                lamMaxm4 = DlamMax1Dx * (unm4 + cm4);
                lamMaxm5 = DlamMax1Dx * (unm5 + cm5);
            }
        }
        else {
            lamMax = - (un-Vgn) + c - minTau;
            DlamMax1Dx = 2*(atan(b*lamMax)/pi + 1.0/2.0 + (b*lamMax) / ((1+b*lamMax*b*lamMax) * pi)) - 1.0;
            lamMax = 2*(lamMax*(atan(b*lamMax)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0) - lamMax + minTau;
            if (computeJacobian == 1) {
                lamMaxm1 = DlamMax1Dx * (-unm1 + cm1);
                lamMaxm2 = DlamMax1Dx * (-unm2 + cm2);
                lamMaxm3 = DlamMax1Dx * (-unm3 + cm3);
                lamMaxm4 = DlamMax1Dx * (-unm4 + cm4);
                lamMaxm5 = DlamMax1Dx * (-unm5 + cm5);
            }
        }
        
       if (epslm > 0) {
           lamMax = 0.5 * (lamMax * lamMax / (epslm * c) + epslm * c);
           lamMaxm1 = lamMax * lamMaxm1 / (epslm*c) + 0.5 * (1 - lamMax * lamMax / (epslm * c * epslm * c)) * (epslm * cm1);
           lamMaxm2 = lamMax * lamMaxm2 / (epslm*c) + 0.5 * (1 - lamMax * lamMax / (epslm * c * epslm * c)) * (epslm * cm2);
           lamMaxm3 = lamMax * lamMaxm3 / (epslm*c) + 0.5 * (1 - lamMax * lamMax / (epslm * c * epslm * c)) * (epslm * cm3);
           lamMaxm4 = lamMax * lamMaxm4 / (epslm*c) + 0.5 * (1 - lamMax * lamMax / (epslm * c * epslm * c)) * (epslm * cm4);
           lamMaxm5 = lamMax * lamMaxm5 / (epslm*c) + 0.5 * (1 - lamMax * lamMax / (epslm * c * epslm * c)) * (epslm * cm5);
       }
        
        for (k = 0; k < 5; k++)
            for (j = 0; j < nch; j++)
                tau[i + j * numPoints + k * sz2] = (j == k) ? lamMax : 0.0;

        if (tauValue > 0) {
            stabTurbEq = tauValue;
            if (computeJacobian == 1) {
                for (j = 0; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }
        else {
            stabTurbEq = lamMax;
            if (computeJacobian == 1) {
                DstabTurbEq[0] = lamMaxm1;
                DstabTurbEq[1] = lamMaxm2;
                DstabTurbEq[2] = lamMaxm3;
                DstabTurbEq[3] = lamMaxm4;
                DstabTurbEq[4] = lamMaxm5;
                for (j = 5; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }

        for (k = 5; k < nch; k++)
            for (j = 0; j < nch; j++)
                tau[i + j * numPoints + k * sz2] = (j == k) ? stabTurbEq : 0.0;
        // tau: numPoints / nch / nch

        if (computeJacobian == 1) {
            int sz3 = nch * sz2;
            DlamMax[0] = lamMaxm1;
            DlamMax[1] = lamMaxm2;
            DlamMax[2] = lamMaxm3;
            DlamMax[3] = lamMaxm4;
            DlamMax[4] = lamMaxm5;
            for (j=5; j<nch; j++)
                DlamMax[j] = 0.0;

            for (m = 0; m < nch; m++)
                for (k = 0; k < 5; k++)
                    for (j = 0; j < nch; j++)
                        tau_uh[i + j * numPoints + k * sz2 + m * sz3] = (j == k) ? DlamMax[m] : 0.0;

            for (m = 0; m < nch; m++)
                for (k = 5; k < nch; k++)
                    for (j = 0; j < nch; j++)
                        tau_uh[i + j * numPoints + k * sz2 + m * sz3] = (j == k) ? DstabTurbEq[m] : 0.0;
        }
        // tau_uh: numPoints / nch / nch / nch
    }
    
    for (i=0; i<numPoints*nch*nch*nch; i++)
        tau_uh[i] = 0.0;
}

void getLaxFriedrichStabilizationTensor(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    if (nd == 2)
        getLaxFriedrichStabilizationTensor2d(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (nd == 3)
        getLaxFriedrichStabilizationTensor3d(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else
        error("Number of dimensions not implemented.\n");
    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch
}

void getRoeStabilizationTensor2d(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    // This functions works for compressible Euler, NS and RANS (any turbulence model).
    // If tauValue <= 0, the maximum eigenvalue is used to stabilize the turbulence model equation(s). Otherwise, tauValue is used.

//    getan2d(tau, tau_uh, uh, nl, param, 1, numPoints);

    double gam, gam1, epslm, Minf, tauValue, signun, ic, Dc2Reg_Dc2;
    double nx, ny;
    double r, ru, rv, rE;
    double r1, r1m1, r1m2, r1m3, r1m4;
    double uv, uvm1, uvm2, uvm3, uvm4;
    double vv, vvm1, vvm2, vvm3, vvm4;
    double Vgx, Vgy, Vgn;
    double E, Em1, Em2, Em3, Em4;
    double af, afm1, afm2, afm3, afm4;
    double p, pm1, pm2, pm3, pm4;
    double h, hm1, hm2, hm3, hm4;
    double s1, s1m1, s1m2, s1m3, s1m4;
    double s2, s2m1, s2m2, s2m3, s2m4;
    double c2, c2m1, c2m2, c2m3, c2m4;
    double c, cm1, cm2, cm3, cm4;
    double un, unm1, unm2, unm3, unm4;
    double cc1, cc1m1, cc1m2, cc1m3, cc1m4;
    double cc2, cc2m1, cc2m2, cc2m3, cc2m4;
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4;
    double rlam1, rlam1m1, rlam1m2, rlam1m3, rlam1m4;
    double rlam2, rlam2m1, rlam2m2, rlam2m3, rlam2m4;
    double rlam3, rlam3m1, rlam3m2, rlam3m3, rlam3m4;
    double lamMax, lamMaxm1, lamMaxm2, lamMaxm3, lamMaxm4;
    
    gam   = param[0];
    gam1 = gam - 1.0;
    epslm = param[1];
    Minf   = param[4];
    tauValue = param[5];
    
    double x, Drlam1Dx, Drlam2Dx, Drlam3Dx;
    double pi = 3.141592653589793;
    double b = 100.0;
    double minTau = 0.1;
    double rMin = 0.5;
    double pMin = (1/(gam*Minf*Minf)) / 2.0;//1.0e-3;
    double cMin = 1.0e-3;
    
    double stabTurbEq;
    double *DstabTurbEq  = new double [nch];
    
    int nch2 = nch * nch, sz2 = nch * numPoints, sz3 = sz2 * nch;
    int i, j, k, m, index_tmp;

    Int ALEflag = app.ALEflag;

    for (i=0; i<numPoints; i++) {
        Int useConstantTau = 0;
        
        r   = uh[0*numPoints+i]; // 1.0000;
        ru  = uh[1*numPoints+i]; // 0.9994;
        rv  = uh[2*numPoints+i]; // 0.0349;
        rE  = uh[3*numPoints+i]; // 45.1429;
        
        ////////////////////////////////
        ////////////////////////////////
        // Fix to improve robustness near non-physical states:
//         if (r < rMin) {
//             rE = abs(rMin * (rE / r));
//             r = rMin;
//         }
//         r1   = 1.0/r;
//         uv   = ru*r1;
//         vv   = rv*r1;
//         E    = rE*r1;
//         af   = 0.5*(uv*uv+vv*vv);
//         p    = gam1*(rE -r*af);
//         if (p < pMin) {
//             rE = pMin / gam1 + r*af;
//         }
        ////////////////////////////////
        ////////////////////////////////
        
        nx  = nl[0*numPoints+i];
        ny  = nl[1*numPoints+i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
        }

        r1   = 1.0/r;
        uv   = ru*r1;
        vv   = rv*r1;
        E    = rE*r1;
        af   = 0.5*(uv*uv+vv*vv);
        p    = gam1*(rE -r*af);
        h    = E   + p*r1;
        c2   = gam* p*r1;
        
        if (r < rMin || p < pMin || c2 < cMin*cMin)
            useConstantTau = 1;
        
        Dc2Reg_Dc2 = atan(b*c2)/pi + 1.0/2.0 + (b*c2) / ((1+b*c2*b*c2) * pi);
        c2 = max ( c2*(atan(b*c2)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 , cMin*cMin);
        c    = sqrt(c2);
        un   = uv*nx   + vv*ny;

        r1m1 = -1.0/(r*r);
        r1m2 = 0.0;
        r1m3 = 0.0;
        r1m4 = 0.0;

        if (computeJacobian == 1) {
            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;
            uvm3 = ru * r1m3;
            uvm4 = ru * r1m4;

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;
            vvm3 = r1 + rv * r1m3;
            vvm4 = rv * r1m4;

            Em1 = rE * r1m1;
            Em2 = rE * r1m2;
            Em3 = rE * r1m3;
            Em4 = r1 + rE * r1m4;

            afm1 = uv * uvm1 + vv * vvm1;
            afm2 = uv * uvm2 + vv * vvm2;
            afm3 = uv * uvm3 + vv * vvm3;
            afm4 = uv * uvm4 + vv * vvm4;

            pm1 = gam1 * (-af - r * afm1);
            pm2 = gam1 * (-r * afm2);
            pm3 = gam1 * (-r * afm3);
            pm4 = gam1 * (1.0 - r * afm4);

            hm1 = Em1 + pm1 * r1 + p * r1m1;
            hm2 = Em2 + pm2 * r1 + p * r1m2;
            hm3 = Em3 + pm3 * r1 + p * r1m3;
            hm4 = Em4 + pm4 * r1 + p * r1m4;

            c2m1 = gam * (pm1 * r1 + p * r1m1);
            c2m2 = gam * (pm2 * r1 + p * r1m2);
            c2m3 = gam * (pm3 * r1 + p * r1m3);
            c2m4 = gam * (pm4 * r1 + p * r1m4);
            
            c2m1 *= Dc2Reg_Dc2;
            c2m2 *= Dc2Reg_Dc2;
            c2m3 *= Dc2Reg_Dc2;
            c2m4 *= Dc2Reg_Dc2;
            
            cm1 = 0.5 * c2m1 / c;
            cm2 = 0.5 * c2m2 / c;
            cm3 = 0.5 * c2m3 / c;
            cm4 = 0.5 * c2m4 / c;

            unm1 = uvm1 * nx + vvm1 * ny;
            unm2 = uvm2 * nx + vvm2 * ny;
            unm3 = uvm3 * nx + vvm3 * ny;
            unm4 = uvm4 * nx + vvm4 * ny;
        }

        Vgn = Vgx * nx + Vgy * ny;
        
        x = un - Vgn + c - minTau;
        rlam1 = 2*(x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0) - x + minTau;
        Drlam1Dx = 2*(atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi)) - 1.0;
        if (computeJacobian == 1) {
            rlam1m1 = Drlam1Dx * (unm1 + cm1);
            rlam1m2 = Drlam1Dx * (unm2 + cm2);
            rlam1m3 = Drlam1Dx * (unm3 + cm3);
            rlam1m4 = Drlam1Dx * (unm4 + cm4);
        }

        if (epslm > 0) {
            rlam = 0.5 * (rlam1 * rlam1 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam1 * rlam1m1 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam1 * rlam1m2 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam1 * rlam1m3 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam1 * rlam1m4 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm4);
            }

            ic = rlam1 < epslm * c;
            rlam1 = ic * rlam + (1 - ic) * rlam1;

            if (computeJacobian == 1) {
                rlam1m1 = ic * rlamm1 + (1 - ic) * rlam1m1;
                rlam1m2 = ic * rlamm2 + (1 - ic) * rlam1m2;
                rlam1m3 = ic * rlamm3 + (1 - ic) * rlam1m3;
                rlam1m4 = ic * rlamm4 + (1 - ic) * rlam1m4;
            }
        }

        x = un - Vgn - c - minTau;
        rlam2 = 2*(x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0) - x + minTau;
        Drlam2Dx = 2*(atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi)) - 1.0;
        if (computeJacobian == 1) {
            rlam2m1 = Drlam2Dx * (unm1 - cm1);
            rlam2m2 = Drlam2Dx * (unm2 - cm2);
            rlam2m3 = Drlam2Dx * (unm3 - cm3);
            rlam2m4 = Drlam2Dx * (unm4 - cm4);
        }

        if (epslm > 0) {
            rlam = 0.5 * (rlam2 * rlam2 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam2 * rlam2m1 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam2 * rlam2m2 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam2 * rlam2m3 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam2 * rlam2m4 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm4);
            }

            ic = rlam2 < epslm * c;
            rlam2 = ic * rlam + (1 - ic) * rlam2;

            if (computeJacobian == 1) {
                rlam2m1 = ic * rlamm1 + (1 - ic) * rlam2m1;
                rlam2m2 = ic * rlamm2 + (1 - ic) * rlam2m2;
                rlam2m3 = ic * rlamm3 + (1 - ic) * rlam2m3;
                rlam2m4 = ic * rlamm4 + (1 - ic) * rlam2m4;
            }
        }

        x = un - Vgn - minTau;
        rlam3 = 2*(x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0) - x + minTau;
        Drlam3Dx = 2*(atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi)) - 1.0;
        if (computeJacobian == 1) {
            rlam3m1 = Drlam3Dx * unm1;
            rlam3m2 = Drlam3Dx * unm2;
            rlam3m3 = Drlam3Dx * unm3;
            rlam3m4 = Drlam3Dx * unm4;
        }

        if (epslm > 0) {
            rlam = 0.5 * (rlam3 * rlam3 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam3 * rlam3m1 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam3 * rlam3m2 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam3 * rlam3m3 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam3 * rlam3m4 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm4);
            }

            ic = rlam3 < epslm * c;
            rlam3 = ic * rlam + (1 - ic) * rlam3;

            if (computeJacobian == 1) {
                rlam3m1 = ic * rlamm1 + (1 - ic) * rlam3m1;
                rlam3m2 = ic * rlamm2 + (1 - ic) * rlam3m2;
                rlam3m3 = ic * rlamm3 + (1 - ic) * rlam3m3;
                rlam3m4 = ic * rlamm4 + (1 - ic) * rlam3m4;
            }
        }

        if (un - Vgn >= 0) {
            lamMax = (un-Vgn) + c;
            if (computeJacobian == 1) {
                lamMaxm1 = unm1 + cm1;
                lamMaxm2 = unm2 + cm2;
                lamMaxm3 = unm3 + cm3;
                lamMaxm4 = unm4 + cm4;
            }
        }
        else {
            lamMax = - (un-Vgn) + c;
            if (computeJacobian == 1) {
                lamMaxm1 = -unm1 + cm1;
                lamMaxm2 = -unm2 + cm2;
                lamMaxm3 = -unm3 + cm3;
                lamMaxm4 = -unm4 + cm4;
            }
        }

        if (tauValue > 0) {
            stabTurbEq = tauValue;
            if (computeJacobian == 1) {
                for (j = 0; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }
        else {
            stabTurbEq = lamMax;
            if (computeJacobian == 1) {
                DstabTurbEq[0] = lamMaxm1;
                DstabTurbEq[1] = lamMaxm2;
                DstabTurbEq[2] = lamMaxm3;
                DstabTurbEq[3] = lamMaxm4;
                for (j = 4; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }

        s1      = 0.5*(rlam1   + rlam2);
        s2      = 0.5*(rlam1   - rlam2);

        if (computeJacobian == 1) {
            s1m1 = 0.5 * (rlam1m1 + rlam2m1);
            s1m2 = 0.5 * (rlam1m2 + rlam2m2);
            s1m3 = 0.5 * (rlam1m3 + rlam2m3);
            s1m4 = 0.5 * (rlam1m4 + rlam2m4);

            s2m1 = 0.5 * (rlam1m1 - rlam2m1);
            s2m2 = 0.5 * (rlam1m2 - rlam2m2);
            s2m3 = 0.5 * (rlam1m3 - rlam2m3);
            s2m4 = 0.5 * (rlam1m4 - rlam2m4);
        }

        cc1   = gam1*(s1-rlam3)*af/c2-(s2*un/c);
        cc2   = gam1*s2*af/c-(s1-rlam3)*un;

        index_tmp = 0*nch*numPoints+i;
        tau[index_tmp + 0*numPoints]  = rlam3+cc1;
        tau[index_tmp + 1*numPoints]  = cc1*uv+cc2*nx;
        tau[index_tmp + 2*numPoints]  = cc1*vv+cc2*ny;
        tau[index_tmp + 3*numPoints]  = cc1*h+cc2*un;
        for (j=4; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = gam1*((s1m1-rlam3m1)*af/c2 + (s1-rlam3)*afm1/c2 - (s1-rlam3)*af*c2m1/(c2*c2)) - s2m1*un/c - s2*unm1/c + s2*un*cm1/(c*c);
            cc1m2 = gam1*((s1m2-rlam3m2)*af/c2 + (s1-rlam3)*afm2/c2 - (s1-rlam3)*af*c2m2/(c2*c2)) - s2m2*un/c - s2*unm2/c + s2*un*cm2/(c*c);
            cc1m3 = gam1*((s1m3-rlam3m3)*af/c2 + (s1-rlam3)*afm3/c2 - (s1-rlam3)*af*c2m3/(c2*c2)) - s2m3*un/c - s2*unm3/c + s2*un*cm3/(c*c);
            cc1m4 = gam1*((s1m4-rlam3m4)*af/c2 + (s1-rlam3)*afm4/c2 - (s1-rlam3)*af*c2m4/(c2*c2)) - s2m4*un/c - s2*unm4/c + s2*un*cm4/(c*c);

            cc2m1 = gam1*(s2m1*af/c + s2*afm1/c - s2*af*cm1/(c*c)) - (s1m1-rlam3m1)*un - (s1-rlam3)*unm1;
            cc2m2 = gam1*(s2m2*af/c + s2*afm2/c - s2*af*cm2/(c*c)) - (s1m2-rlam3m2)*un - (s1-rlam3)*unm2;
            cc2m3 = gam1*(s2m3*af/c + s2*afm3/c - s2*af*cm3/(c*c)) - (s1m3-rlam3m3)*un - (s1-rlam3)*unm3;
            cc2m4 = gam1*(s2m4*af/c + s2*afm4/c - s2*af*cm4/(c*c)) - (s1m4-rlam3m4)*un - (s1-rlam3)*unm4;

            index_tmp = (0*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m1 + cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m2 + cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m3 + cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m4 + cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;

            index_tmp = 0*sz2 + i;
            for (k=4; k<nch; k++)
                for (j=0; j<4; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=4; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1   = -gam1*(s1-rlam3)*uv/c2+(s2*nx/c);
        cc2   = -gam1*s2*uv/c + (s1-rlam3)*nx;

        index_tmp = 1*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = rlam3+cc1*uv+cc2*nx;
        tau[index_tmp + 2*numPoints] = cc1*vv+cc2*ny;
        tau[index_tmp + 3*numPoints] = cc1*h+cc2*un;
        for (j=4; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = -gam1*((s1m1-rlam3m1)*uv/c2 + (s1-rlam3)*uvm1/c2 - (s1-rlam3)*uv*c2m1/(c2*c2)) + s2m1*nx/c - s2*nx*cm1/(c*c);
            cc1m2 = -gam1*((s1m2-rlam3m2)*uv/c2 + (s1-rlam3)*uvm2/c2 - (s1-rlam3)*uv*c2m2/(c2*c2)) + s2m2*nx/c - s2*nx*cm2/(c*c);
            cc1m3 = -gam1*((s1m3-rlam3m3)*uv/c2 + (s1-rlam3)*uvm3/c2 - (s1-rlam3)*uv*c2m3/(c2*c2)) + s2m3*nx/c - s2*nx*cm3/(c*c);
            cc1m4 = -gam1*((s1m4-rlam3m4)*uv/c2 + (s1-rlam3)*uvm4/c2 - (s1-rlam3)*uv*c2m4/(c2*c2)) + s2m4*nx/c - s2*nx*cm4/(c*c);

            cc2m1 = -gam1*(s2m1*uv/c + s2*uvm1/c - s2*uv*cm1/(c*c)) + (s1m1-rlam3m1)*nx;
            cc2m2 = -gam1*(s2m2*uv/c + s2*uvm2/c - s2*uv*cm2/(c*c)) + (s1m2-rlam3m2)*nx;
            cc2m3 = -gam1*(s2m3*uv/c + s2*uvm3/c - s2*uv*cm3/(c*c)) + (s1m3-rlam3m3)*nx;
            cc2m4 = -gam1*(s2m4*uv/c + s2*uvm4/c - s2*uv*cm4/(c*c)) + (s1m4-rlam3m4)*nx;

            index_tmp = (0*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = rlam3m1 + cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = rlam3m2 + cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = rlam3m3 + cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = rlam3m4 + cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;

            index_tmp = 1*sz2 + i;
            for (k=4; k<nch; k++)
                for (j=0; j<4; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=4; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1   = -gam1*(s1-rlam3)*vv/c2+(s2*ny/c);
        cc2   = -gam1*s2*vv/c+(s1-rlam3)*ny;

        index_tmp = 2*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = cc1*uv+cc2*nx;
        tau[index_tmp + 2*numPoints] = rlam3+cc1*vv+cc2*ny;
        tau[index_tmp + 3*numPoints] = cc1*h+cc2*un;
        for (j=4; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = -gam1*((s1m1-rlam3m1)*vv/c2 + (s1-rlam3)*vvm1/c2 - (s1-rlam3)*vv*c2m1/(c2*c2)) + s2m1*ny/c - s2*ny*cm1/(c*c);
            cc1m2 = -gam1*((s1m2-rlam3m2)*vv/c2 + (s1-rlam3)*vvm2/c2 - (s1-rlam3)*vv*c2m2/(c2*c2)) + s2m2*ny/c - s2*ny*cm2/(c*c);
            cc1m3 = -gam1*((s1m3-rlam3m3)*vv/c2 + (s1-rlam3)*vvm3/c2 - (s1-rlam3)*vv*c2m3/(c2*c2)) + s2m3*ny/c - s2*ny*cm3/(c*c);
            cc1m4 = -gam1*((s1m4-rlam3m4)*vv/c2 + (s1-rlam3)*vvm4/c2 - (s1-rlam3)*vv*c2m4/(c2*c2)) + s2m4*ny/c - s2*ny*cm4/(c*c);

            cc2m1 = -gam1*(s2m1*vv/c + s2*vvm1/c - s2*vv*cm1/(c*c)) + (s1m1-rlam3m1)*ny;
            cc2m2 = -gam1*(s2m2*vv/c + s2*vvm2/c - s2*vv*cm2/(c*c)) + (s1m2-rlam3m2)*ny;
            cc2m3 = -gam1*(s2m3*vv/c + s2*vvm3/c - s2*vv*cm3/(c*c)) + (s1m3-rlam3m3)*ny;
            cc2m4 = -gam1*(s2m4*vv/c + s2*vvm4/c - s2*vv*cm4/(c*c)) + (s1m4-rlam3m4)*ny;

            index_tmp = (0*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m1 + cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m2 + cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m3 + cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m4 + cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;

            index_tmp = 2*sz2 + i;
            for (k=4; k<nch; k++)
                for (j=0; j<4; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=4; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1   = gam1*(s1-rlam3)/c2;
        cc2   = gam1*s2/c;

        index_tmp = 3*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = cc1*uv+cc2*nx;
        tau[index_tmp + 2*numPoints] = cc1*vv+cc2*ny;
        tau[index_tmp + 3*numPoints] = rlam3+cc1*h+cc2*un;
        for (j=4; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = gam1*((s1m1-rlam3m1)/c2 - (s1-rlam3)*c2m1/(c2*c2));
            cc1m2 = gam1*((s1m2-rlam3m2)/c2 - (s1-rlam3)*c2m2/(c2*c2));
            cc1m3 = gam1*((s1m3-rlam3m3)/c2 - (s1-rlam3)*c2m3/(c2*c2));
            cc1m4 = gam1*((s1m4-rlam3m4)/c2 - (s1-rlam3)*c2m4/(c2*c2));

            cc2m1 = gam1*(s2m1/c - s2*cm1/(c*c));
            cc2m2 = gam1*(s2m2/c - s2*cm2/(c*c));
            cc2m3 = gam1*(s2m3/c - s2*cm3/(c*c));
            cc2m4 = gam1*(s2m4/c - s2*cm4/(c*c));

            index_tmp = (0*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m1 + cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m2 + cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m3 + cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m4 + cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;

            index_tmp = 3*sz2 + i;
            for (k=4; k<nch; k++)
                for (j=0; j<4; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=4; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        for (k=4; k<nch; k++)
            for (j=0; j<nch; j++)
                tau[k*sz2 + j*numPoints + i] = (k == j) ? stabTurbEq : 0.0;

        if (computeJacobian == 1) {
            for (m = 0; m < nch; m++)
                for (k = 4; k < nch; k++)
                    for (j = 0; j < nch; j++)
                        tau_uh[m*sz3 + k*sz2 + j*numPoints + i] = (k == j) ? DstabTurbEq[m] : 0.0;
        }
        
        if (useConstantTau == 1) {         
            for (k=0; k<nch; k++)
                for (j=0; j<nch; j++)
                    for (i = 0; i < numPoints; i++)
                        tau[i + j * numPoints + k * sz2] = (j==k) ? tauValue : 0.0;
                        
            if (computeJacobian == 1) {
                for (i=0; i<numPoints * nch * nch * nch; i++)
                    tau_uh[i] = 0.0;
            }    
        }
    }
    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch
    
   for (i=0; i<numPoints*nch*nch*nch; i++)
       tau_uh[i] = 0.0;

    delete[] DstabTurbEq;
}

void getRoeStabilizationTensor3d(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    // This functions works for compressible Euler, NS and RANS (any turbulence model).
    // If tauValue <= 0, the maximum eigenvalue is used to stabilize the turbulence model equation(s). Otherwise, tauValue is used.

//    getan3d(tau, tau_uh, uh, nl, param, 0, numPoints);

    double gam, gam1, epslm, tauValue, signun, ic, Dc2Reg_Dc2;
    double nx, ny, nz;
    double r, ru, rv, rw, rE;
    double r1, r1m1, r1m2, r1m3, r1m4, r1m5;
    double uv, uvm1, uvm2, uvm3, uvm4, uvm5;
    double vv, vvm1, vvm2, vvm3, vvm4, vvm5;
    double wv, wvm1, wvm2, wvm3, wvm4, wvm5;
    double Vgx, Vgy, Vgz, Vgn;
    double E, Em1, Em2, Em3, Em4, Em5;
    double af, afm1, afm2, afm3, afm4, afm5;
    double p, pm1, pm2, pm3, pm4, pm5;
    double h, hm1, hm2, hm3, hm4, hm5;
    double s1, s1m1, s1m2, s1m3, s1m4, s1m5;
    double s2, s2m1, s2m2, s2m3, s2m4, s2m5;
    double c2, c2m1, c2m2, c2m3, c2m4, c2m5;
    double c, cm1, cm2, cm3, cm4, cm5;
    double un, unm1, unm2, unm3, unm4, unm5;
    double cc1, cc1m1, cc1m2, cc1m3, cc1m4, cc1m5;
    double cc2, cc2m1, cc2m2, cc2m3, cc2m4, cc2m5;
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4, rlamm5;
    double rlam1, rlam1m1, rlam1m2, rlam1m3, rlam1m4, rlam1m5;
    double rlam2, rlam2m1, rlam2m2, rlam2m3, rlam2m4, rlam2m5;
    double rlam3, rlam3m1, rlam3m2, rlam3m3, rlam3m4, rlam3m5;
    double lamMax, lamMaxm1, lamMaxm2, lamMaxm3, lamMaxm4, lamMaxm5;

    gam = param[0];
    gam1 = gam - 1.0;
    epslm = param[1];
    tauValue = param[5];

    double x, Drlam1Dx, Drlam2Dx, Drlam3Dx;
    double pi = 3.141592653589793;
    double b = 100.0;
    double minTau = 0.1;
    
    double stabTurbEq;
    double *DstabTurbEq  = new double [nch];

    int nch2 = nch * nch, sz2 = nch * numPoints, sz3 = sz2 * nch;
    int i, j, k, m, index_tmp;

    Int ALEflag = app.ALEflag;

    for (i = 0; i < numPoints; i++) {
        nx = nl[0 * numPoints + i];
        ny = nl[1 * numPoints + i];
        nz = nl[2 * numPoints + i];

        r = uh[0 * numPoints + i];
        ru = uh[1 * numPoints + i];
        rv = uh[2 * numPoints + i];
        rw = uh[3 * numPoints + i];
        rE = uh[4 * numPoints + i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
            Vgz = pg[(2 * nd + 2) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
            Vgz = 0.0;
        }

        r1 = 1.0 / r;
        uv = ru * r1;
        vv = rv * r1;
        wv = rw * r1;
        E = rE * r1;
        af = 0.5 * (uv * uv + vv * vv + wv * wv);
        p = gam1 * (rE - r * af);
        h = E + p * r1;
        c2 = gam * p * r1;
        Dc2Reg_Dc2 = atan(b*c2)/pi + 1.0/2.0 + (b*c2) / ((1+b*c2*b*c2) * pi);
        c2 = max ( c2*(atan(b*c2)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 , 1.0e-6);
        c = sqrt(c2);
        un = uv * nx + vv * ny + wv * nz;

        if (computeJacobian == 1) {
            r1m1 = -1.0 / (r * r);
            r1m2 = 0.0;
            r1m3 = 0.0;
            r1m4 = 0.0;
            r1m5 = 0.0;

            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;
            uvm3 = ru * r1m3;
            uvm4 = ru * r1m4;
            uvm5 = ru * r1m5;

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;
            vvm3 = r1 + rv * r1m3;
            vvm4 = rv * r1m4;
            vvm5 = rv * r1m5;

            wvm1 = rw * r1m1;
            wvm2 = rw * r1m2;
            wvm3 = rw * r1m3;
            wvm4 = r1 + rw * r1m4;
            wvm5 = rw * r1m5;

            Em1 = rE * r1m1;
            Em2 = rE * r1m2;
            Em3 = rE * r1m3;
            Em4 = rE * r1m4;
            Em5 = r1 + rE * r1m5;

            afm1 = uv * uvm1 + vv * vvm1 + wv * wvm1;
            afm2 = uv * uvm2 + vv * vvm2 + wv * wvm2;
            afm3 = uv * uvm3 + vv * vvm3 + wv * wvm3;
            afm4 = uv * uvm4 + vv * vvm4 + wv * wvm4;
            afm5 = uv * uvm5 + vv * vvm5 + wv * wvm5;

            pm1 = gam1 * (-af - r * afm1);
            pm2 = gam1 * (-r * afm2);
            pm3 = gam1 * (-r * afm3);
            pm4 = gam1 * (-r * afm4);
            pm5 = gam1 * (1.0 - r * afm5);

            hm1 = Em1 + pm1 * r1 + p * r1m1;
            hm2 = Em2 + pm2 * r1 + p * r1m2;
            hm3 = Em3 + pm3 * r1 + p * r1m3;
            hm4 = Em4 + pm4 * r1 + p * r1m4;
            hm5 = Em5 + pm5 * r1 + p * r1m5;

            c2m1 = gam * (pm1 * r1 + p * r1m1);
            c2m2 = gam * (pm2 * r1 + p * r1m2);
            c2m3 = gam * (pm3 * r1 + p * r1m3);
            c2m4 = gam * (pm4 * r1 + p * r1m4);
            c2m5 = gam * (pm5 * r1 + p * r1m5);
            
            c2m1 *= Dc2Reg_Dc2;
            c2m2 *= Dc2Reg_Dc2;
            c2m3 *= Dc2Reg_Dc2;
            c2m4 *= Dc2Reg_Dc2;
            c2m5 *= Dc2Reg_Dc2;

            cm1 = 0.5 * c2m1 / c;
            cm2 = 0.5 * c2m2 / c;
            cm3 = 0.5 * c2m3 / c;
            cm4 = 0.5 * c2m4 / c;
            cm5 = 0.5 * c2m5 / c;

            unm1 = uvm1 * nx + vvm1 * ny + wvm1 * nz;
            unm2 = uvm2 * nx + vvm2 * ny + wvm2 * nz;
            unm3 = uvm3 * nx + vvm3 * ny + wvm3 * nz;
            unm4 = uvm4 * nx + vvm4 * ny + wvm4 * nz;
            unm5 = uvm5 * nx + vvm5 * ny + wvm5 * nz;
        }

        Vgn = Vgx * nx + Vgy * ny + Vgz * nz;
        
        x = un - Vgn + c - minTau;
        rlam1 = 2*(x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0) - x + minTau;
        Drlam1Dx = 2*(atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi)) - 1.0;
        if (computeJacobian == 1) {
            rlam1m1 = Drlam1Dx * (unm1 + cm1);
            rlam1m2 = Drlam1Dx * (unm2 + cm2);
            rlam1m3 = Drlam1Dx * (unm3 + cm3);
            rlam1m4 = Drlam1Dx * (unm4 + cm4);
            rlam1m5 = Drlam1Dx * (unm5 + cm5);
        }

        if (epslm > 0) {
            rlam = 0.5 * (rlam1 * rlam1 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam1 * rlam1m1 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam1 * rlam1m2 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam1 * rlam1m3 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam1 * rlam1m4 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm4);
                rlamm5 =
                        rlam1 * rlam1m5 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm5);
            }

            ic = rlam1 < epslm * c;
            rlam1 = ic * rlam + (1 - ic) * rlam1;

            if (computeJacobian == 1) {
                rlam1m1 = ic * rlamm1 + (1 - ic) * rlam1m1;
                rlam1m2 = ic * rlamm2 + (1 - ic) * rlam1m2;
                rlam1m3 = ic * rlamm3 + (1 - ic) * rlam1m3;
                rlam1m4 = ic * rlamm4 + (1 - ic) * rlam1m4;
                rlam1m5 = ic * rlamm5 + (1 - ic) * rlam1m5;
            }
        }
        
        x = un - Vgn - c - minTau;
        rlam2 = 2*(x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0) - x + minTau;
        Drlam2Dx = 2*(atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi)) - 1.0;
        if (computeJacobian == 1) {
            rlam2m1 = Drlam2Dx * (unm1 - cm1);
            rlam2m2 = Drlam2Dx * (unm2 - cm2);
            rlam2m3 = Drlam2Dx * (unm3 - cm3);
            rlam2m4 = Drlam2Dx * (unm4 - cm4);
            rlam2m5 = Drlam2Dx * (unm5 - cm5);
        }
        
        if (epslm > 0) {
            rlam = 0.5 * (rlam2 * rlam2 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam2 * rlam2m1 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam2 * rlam2m2 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam2 * rlam2m3 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam2 * rlam2m4 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm4);
                rlamm5 =
                        rlam2 * rlam2m5 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm5);
            }

            ic = rlam2 < epslm * c;
            rlam2 = ic * rlam + (1 - ic) * rlam2;

            if (computeJacobian == 1) {
                rlam2m1 = ic * rlamm1 + (1 - ic) * rlam2m1;
                rlam2m2 = ic * rlamm2 + (1 - ic) * rlam2m2;
                rlam2m3 = ic * rlamm3 + (1 - ic) * rlam2m3;
                rlam2m4 = ic * rlamm4 + (1 - ic) * rlam2m4;
                rlam2m5 = ic * rlamm5 + (1 - ic) * rlam2m5;
            }
        }
        
        x = un - Vgn - minTau;
        rlam3 = 2*(x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0) - x + minTau;
        Drlam3Dx = 2*(atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi)) - 1.0;
        if (computeJacobian == 1) {
            rlam3m1 = Drlam3Dx * unm1;
            rlam3m2 = Drlam3Dx * unm2;
            rlam3m3 = Drlam3Dx * unm3;
            rlam3m4 = Drlam3Dx * unm4;
            rlam3m5 = Drlam3Dx * unm5;
        }

        if (epslm > 0) {
            rlam = 0.5 * (rlam3 * rlam3 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam3 * rlam3m1 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam3 * rlam3m2 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam3 * rlam3m3 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam3 * rlam3m4 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm4);
                rlamm5 =
                        rlam3 * rlam3m5 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm5);
            }

            ic = rlam3 < epslm * c;
            rlam3 = ic * rlam + (1 - ic) * rlam3;

            if (computeJacobian == 1) {
                rlam3m1 = ic * rlamm1 + (1 - ic) * rlam3m1;
                rlam3m2 = ic * rlamm2 + (1 - ic) * rlam3m2;
                rlam3m3 = ic * rlamm3 + (1 - ic) * rlam3m3;
                rlam3m4 = ic * rlamm4 + (1 - ic) * rlam3m4;
                rlam3m5 = ic * rlamm5 + (1 - ic) * rlam3m5;
            }
        }

        if (un - Vgn >= 0) {
            lamMax = (un-Vgn) + c;
            if (computeJacobian == 1) {
                lamMaxm1 = unm1 + cm1;
                lamMaxm2 = unm2 + cm2;
                lamMaxm3 = unm3 + cm3;
                lamMaxm4 = unm4 + cm4;
                lamMaxm5 = unm5 + cm5;
            }
        }
        else {
            lamMax = - (un-Vgn) + c;
            if (computeJacobian == 1) {
                lamMaxm1 = -unm1 + cm1;
                lamMaxm2 = -unm2 + cm2;
                lamMaxm3 = -unm3 + cm3;
                lamMaxm4 = -unm4 + cm4;
                lamMaxm5 = -unm5 + cm5;
            }
        }

        if (tauValue > 0) {
            stabTurbEq = tauValue;
            if (computeJacobian == 1) {
                for (j = 0; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }
        else {
            stabTurbEq = lamMax;
            if (computeJacobian == 1) {
                DstabTurbEq[0] = lamMaxm1;
                DstabTurbEq[1] = lamMaxm2;
                DstabTurbEq[2] = lamMaxm3;
                DstabTurbEq[3] = lamMaxm4;
                DstabTurbEq[4] = lamMaxm5;
                for (j = 5; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }

        s1 = 0.5 * (rlam1 + rlam2);
        s2 = 0.5 * (rlam1 - rlam2);

        if (computeJacobian == 1) {
            s1m1 = 0.5 * (rlam1m1 + rlam2m1);
            s1m2 = 0.5 * (rlam1m2 + rlam2m2);
            s1m3 = 0.5 * (rlam1m3 + rlam2m3);
            s1m4 = 0.5 * (rlam1m4 + rlam2m4);
            s1m5 = 0.5 * (rlam1m5 + rlam2m5);

            s2m1 = 0.5 * (rlam1m1 - rlam2m1);
            s2m2 = 0.5 * (rlam1m2 - rlam2m2);
            s2m3 = 0.5 * (rlam1m3 - rlam2m3);
            s2m4 = 0.5 * (rlam1m4 - rlam2m4);
            s2m5 = 0.5 * (rlam1m5 - rlam2m5);
        }

        cc1 = gam1 * (s1 - rlam3) * af / c2 - (s2 * un / c);
        cc2 = gam1 * s2 * af / c - (s1 - rlam3) * un;

        index_tmp = 0*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = rlam3 + cc1;
        tau[index_tmp + 1*numPoints] = cc1 * uv + cc2 * nx;
        tau[index_tmp + 2*numPoints] = cc1 * vv + cc2 * ny;
        tau[index_tmp + 3*numPoints] = cc1 * wv + cc2 * nz;
        tau[index_tmp + 4*numPoints] = cc1 * h + cc2 * un;
        for (j=5; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = gam1 * ((s1m1 - rlam3m1) * af / c2 + (s1 - rlam3) * afm1 / c2 - (s1 - rlam3) * af * c2m1 / (c2 * c2)) -
                    s2m1 * un / c - s2 * unm1 / c + s2 * un * cm1 / (c * c);
            cc1m2 = gam1 * ((s1m2 - rlam3m2) * af / c2 + (s1 - rlam3) * afm2 / c2 - (s1 - rlam3) * af * c2m2 / (c2 * c2)) -
                    s2m2 * un / c - s2 * unm2 / c + s2 * un * cm2 / (c * c);
            cc1m3 = gam1 * ((s1m3 - rlam3m3) * af / c2 + (s1 - rlam3) * afm3 / c2 - (s1 - rlam3) * af * c2m3 / (c2 * c2)) -
                    s2m3 * un / c - s2 * unm3 / c + s2 * un * cm3 / (c * c);
            cc1m4 = gam1 * ((s1m4 - rlam3m4) * af / c2 + (s1 - rlam3) * afm4 / c2 - (s1 - rlam3) * af * c2m4 / (c2 * c2)) -
                    s2m4 * un / c - s2 * unm4 / c + s2 * un * cm4 / (c * c);
            cc1m5 = gam1 * ((s1m5 - rlam3m5) * af / c2 + (s1 - rlam3) * afm5 / c2 - (s1 - rlam3) * af * c2m5 / (c2 * c2)) -
                    s2m5 * un / c - s2 * unm5 / c + s2 * un * cm5 / (c * c);

            cc2m1 = gam1 * (s2m1 * af / c + s2 * afm1 / c - s2 * af * cm1 / (c * c)) - (s1m1 - rlam3m1) * un -
                    (s1 - rlam3) * unm1;
            cc2m2 = gam1 * (s2m2 * af / c + s2 * afm2 / c - s2 * af * cm2 / (c * c)) - (s1m2 - rlam3m2) * un -
                    (s1 - rlam3) * unm2;
            cc2m3 = gam1 * (s2m3 * af / c + s2 * afm3 / c - s2 * af * cm3 / (c * c)) - (s1m3 - rlam3m3) * un -
                    (s1 - rlam3) * unm3;
            cc2m4 = gam1 * (s2m4 * af / c + s2 * afm4 / c - s2 * af * cm4 / (c * c)) - (s1m4 - rlam3m4) * un -
                    (s1 - rlam3) * unm4;
            cc2m5 = gam1 * (s2m5 * af / c + s2 * afm5 / c - s2 * af * cm5 / (c * c)) - (s1m5 - rlam3m5) * un -
                    (s1 - rlam3) * unm5;

            index_tmp = (0*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m1 + cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m2 + cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m3 + cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m4 + cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            index_tmp += sz3; // (4*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m5 + cc1m5;
            tau_uh[index_tmp + 1*numPoints] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;

            index_tmp = 0*sz2 + i;
            for (k=5; k<nch; k++)
                for (j=0; j<5; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=5; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1 = -gam1 * (s1 - rlam3) * uv / c2 + (s2 * nx / c);
        cc2 = -gam1 * s2 * uv / c + (s1 - rlam3) * nx;

        index_tmp = 1*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = rlam3 + cc1 * uv + cc2 * nx;
        tau[index_tmp + 2*numPoints] = cc1 * vv + cc2 * ny;
        tau[index_tmp + 3*numPoints] = cc1 * wv + cc2 * nz;
        tau[index_tmp + 4*numPoints] = cc1 * h + cc2 * un;
        for (j=5; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = -gam1 * ((s1m1 - rlam3m1) * uv / c2 + (s1 - rlam3) * uvm1 / c2 - (s1 - rlam3) * uv * c2m1 / (c2 * c2)) +
                    s2m1 * nx / c - s2 * nx * cm1 / (c * c);
            cc1m2 = -gam1 * ((s1m2 - rlam3m2) * uv / c2 + (s1 - rlam3) * uvm2 / c2 - (s1 - rlam3) * uv * c2m2 / (c2 * c2)) +
                    s2m2 * nx / c - s2 * nx * cm2 / (c * c);
            cc1m3 = -gam1 * ((s1m3 - rlam3m3) * uv / c2 + (s1 - rlam3) * uvm3 / c2 - (s1 - rlam3) * uv * c2m3 / (c2 * c2)) +
                    s2m3 * nx / c - s2 * nx * cm3 / (c * c);
            cc1m4 = -gam1 * ((s1m4 - rlam3m4) * uv / c2 + (s1 - rlam3) * uvm4 / c2 - (s1 - rlam3) * uv * c2m4 / (c2 * c2)) +
                    s2m4 * nx / c - s2 * nx * cm4 / (c * c);
            cc1m5 = -gam1 * ((s1m5 - rlam3m5) * uv / c2 + (s1 - rlam3) * uvm5 / c2 - (s1 - rlam3) * uv * c2m5 / (c2 * c2)) +
                    s2m5 * nx / c - s2 * nx * cm5 / (c * c);

            cc2m1 = -gam1 * (s2m1 * uv / c + s2 * uvm1 / c - s2 * uv * cm1 / (c * c)) + (s1m1 - rlam3m1) * nx;
            cc2m2 = -gam1 * (s2m2 * uv / c + s2 * uvm2 / c - s2 * uv * cm2 / (c * c)) + (s1m2 - rlam3m2) * nx;
            cc2m3 = -gam1 * (s2m3 * uv / c + s2 * uvm3 / c - s2 * uv * cm3 / (c * c)) + (s1m3 - rlam3m3) * nx;
            cc2m4 = -gam1 * (s2m4 * uv / c + s2 * uvm4 / c - s2 * uv * cm4 / (c * c)) + (s1m4 - rlam3m4) * nx;
            cc2m5 = -gam1 * (s2m5 * uv / c + s2 * uvm5 / c - s2 * uv * cm5 / (c * c)) + (s1m5 - rlam3m5) * nx;

            index_tmp = (0*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = rlam3m1 + cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = rlam3m2 + cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = rlam3m3 + cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = rlam3m4 + cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            index_tmp += sz3; // (4*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m5;
            tau_uh[index_tmp + 1*numPoints] = rlam3m5 + cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;

            index_tmp = 1*sz2 + i;
            for (k=5; k<nch; k++)
                for (j=0; j<5; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=5; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1 = -gam1 * (s1 - rlam3) * vv / c2 + (s2 * ny / c);
        cc2 = -gam1 * s2 * vv / c + (s1 - rlam3) * ny;

        index_tmp = 2*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = cc1 * uv + cc2 * nx;
        tau[index_tmp + 2*numPoints] = rlam3 + cc1 * vv + cc2 * ny;
        tau[index_tmp + 3*numPoints] = cc1 * wv + cc2 * nz;
        tau[index_tmp + 4*numPoints] = cc1 * h + cc2 * un;
        for (j=5; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = -gam1 * ((s1m1 - rlam3m1) * vv / c2 + (s1 - rlam3) * vvm1 / c2 - (s1 - rlam3) * vv * c2m1 / (c2 * c2)) +
                    s2m1 * ny / c - s2 * ny * cm1 / (c * c);
            cc1m2 = -gam1 * ((s1m2 - rlam3m2) * vv / c2 + (s1 - rlam3) * vvm2 / c2 - (s1 - rlam3) * vv * c2m2 / (c2 * c2)) +
                    s2m2 * ny / c - s2 * ny * cm2 / (c * c);
            cc1m3 = -gam1 * ((s1m3 - rlam3m3) * vv / c2 + (s1 - rlam3) * vvm3 / c2 - (s1 - rlam3) * vv * c2m3 / (c2 * c2)) +
                    s2m3 * ny / c - s2 * ny * cm3 / (c * c);
            cc1m4 = -gam1 * ((s1m4 - rlam3m4) * vv / c2 + (s1 - rlam3) * vvm4 / c2 - (s1 - rlam3) * vv * c2m4 / (c2 * c2)) +
                    s2m4 * ny / c - s2 * ny * cm4 / (c * c);
            cc1m5 = -gam1 * ((s1m5 - rlam3m5) * vv / c2 + (s1 - rlam3) * vvm5 / c2 - (s1 - rlam3) * vv * c2m5 / (c2 * c2)) +
                    s2m5 * ny / c - s2 * ny * cm5 / (c * c);

            cc2m1 = -gam1 * (s2m1 * vv / c + s2 * vvm1 / c - s2 * vv * cm1 / (c * c)) + (s1m1 - rlam3m1) * ny;
            cc2m2 = -gam1 * (s2m2 * vv / c + s2 * vvm2 / c - s2 * vv * cm2 / (c * c)) + (s1m2 - rlam3m2) * ny;
            cc2m3 = -gam1 * (s2m3 * vv / c + s2 * vvm3 / c - s2 * vv * cm3 / (c * c)) + (s1m3 - rlam3m3) * ny;
            cc2m4 = -gam1 * (s2m4 * vv / c + s2 * vvm4 / c - s2 * vv * cm4 / (c * c)) + (s1m4 - rlam3m4) * ny;
            cc2m5 = -gam1 * (s2m5 * vv / c + s2 * vvm5 / c - s2 * vv * cm5 / (c * c)) + (s1m5 - rlam3m5) * ny;

            index_tmp = (0*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m1 + cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m2 + cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m3 + cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m4 + cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            index_tmp += sz3; // (4*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m5;
            tau_uh[index_tmp + 1*numPoints] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m5 + cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;

            index_tmp = 2*sz2 + i;
            for (k=5; k<nch; k++)
                for (j=0; j<5; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=5; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1 = -gam1 * (s1 - rlam3) * wv / c2 + (s2 * nz / c);
        cc2 = -gam1 * s2 * wv / c + (s1 - rlam3) * nz;
        
        index_tmp = 3*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = cc1 * uv + cc2 * nx;
        tau[index_tmp + 2*numPoints] = cc1 * vv + cc2 * ny;
        tau[index_tmp + 3*numPoints] = rlam3 + cc1 * wv + cc2 * nz;
        tau[index_tmp + 4*numPoints] = cc1 * h + cc2 * un;
        for (j=5; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = -gam1 * ((s1m1 - rlam3m1) * wv / c2 + (s1 - rlam3) * wvm1 / c2 - (s1 - rlam3) * wv * c2m1 / (c2 * c2)) +
                    s2m1 * nz / c - s2 * nz * cm1 / (c * c);
            cc1m2 = -gam1 * ((s1m2 - rlam3m2) * wv / c2 + (s1 - rlam3) * wvm2 / c2 - (s1 - rlam3) * wv * c2m2 / (c2 * c2)) +
                    s2m2 * nz / c - s2 * nz * cm2 / (c * c);
            cc1m3 = -gam1 * ((s1m3 - rlam3m3) * wv / c2 + (s1 - rlam3) * wvm3 / c2 - (s1 - rlam3) * wv * c2m3 / (c2 * c2)) +
                    s2m3 * nz / c - s2 * nz * cm3 / (c * c);
            cc1m4 = -gam1 * ((s1m4 - rlam3m4) * wv / c2 + (s1 - rlam3) * wvm4 / c2 - (s1 - rlam3) * wv * c2m4 / (c2 * c2)) +
                    s2m4 * nz / c - s2 * nz * cm4 / (c * c);
            cc1m5 = -gam1 * ((s1m5 - rlam3m5) * wv / c2 + (s1 - rlam3) * wvm5 / c2 - (s1 - rlam3) * wv * c2m5 / (c2 * c2)) +
                    s2m5 * nz / c - s2 * nz * cm5 / (c * c);

            cc2m1 = -gam1 * (s2m1 * wv / c + s2 * wvm1 / c - s2 * wv * cm1 / (c * c)) + (s1m1 - rlam3m1) * nz;
            cc2m2 = -gam1 * (s2m2 * wv / c + s2 * wvm2 / c - s2 * wv * cm2 / (c * c)) + (s1m2 - rlam3m2) * nz;
            cc2m3 = -gam1 * (s2m3 * wv / c + s2 * wvm3 / c - s2 * wv * cm3 / (c * c)) + (s1m3 - rlam3m3) * nz;
            cc2m4 = -gam1 * (s2m4 * wv / c + s2 * wvm4 / c - s2 * wv * cm4 / (c * c)) + (s1m4 - rlam3m4) * nz;
            cc2m5 = -gam1 * (s2m5 * wv / c + s2 * wvm5 / c - s2 * wv * cm5 / (c * c)) + (s1m5 - rlam3m5) * nz;

            index_tmp = (0*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m1 + cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m2 + cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m3 + cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m4 + cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            index_tmp += sz3; // (4*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m5;
            tau_uh[index_tmp + 1*numPoints] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m5 + cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;

            index_tmp = 3*sz2 + i;
            for (k=5; k<nch; k++)
                for (j=0; j<5; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=5; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1 = gam1 * (s1 - rlam3) / c2;
        cc2 = gam1 * s2 / c;

        index_tmp = 4*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = cc1 * uv + cc2 * nx;
        tau[index_tmp + 2*numPoints] = cc1 * vv + cc2 * ny;
        tau[index_tmp + 3*numPoints] = cc1 * wv + cc2 * nz;
        tau[index_tmp + 4*numPoints] = rlam3 + cc1 * h + cc2 * un;
        for (j=5; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = gam1 * ((s1m1 - rlam3m1) / c2 - (s1 - rlam3) * c2m1 / (c2 * c2));
            cc1m2 = gam1 * ((s1m2 - rlam3m2) / c2 - (s1 - rlam3) * c2m2 / (c2 * c2));
            cc1m3 = gam1 * ((s1m3 - rlam3m3) / c2 - (s1 - rlam3) * c2m3 / (c2 * c2));
            cc1m4 = gam1 * ((s1m4 - rlam3m4) / c2 - (s1 - rlam3) * c2m4 / (c2 * c2));
            cc1m5 = gam1 * ((s1m5 - rlam3m5) / c2 - (s1 - rlam3) * c2m5 / (c2 * c2));

            cc2m1 = gam1 * (s2m1 / c - s2 * cm1 / (c * c));
            cc2m2 = gam1 * (s2m2 / c - s2 * cm2 / (c * c));
            cc2m3 = gam1 * (s2m3 / c - s2 * cm3 / (c * c));
            cc2m4 = gam1 * (s2m4 / c - s2 * cm4 / (c * c));
            cc2m5 = gam1 * (s2m5 / c - s2 * cm5 / (c * c));

            index_tmp = (0*nch2 + 4*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            tau_uh[index_tmp + 4*numPoints] = rlam3m1 + cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 4*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            tau_uh[index_tmp + 4*numPoints] = rlam3m2 + cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 4*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            tau_uh[index_tmp + 4*numPoints] = rlam3m3 + cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 4*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            tau_uh[index_tmp + 4*numPoints] = rlam3m4 + cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            index_tmp += sz3; // (4*nch2 + 4*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m5;
            tau_uh[index_tmp + 1*numPoints] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            tau_uh[index_tmp + 4*numPoints] = rlam3m5 + cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;

            index_tmp = 4*sz2 + i;
            for (k=5; k<nch; k++)
                for (j=0; j<5; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=5; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        for (k=5; k<nch; k++)
            for (j=0; j<nch; j++)
                tau[k*sz2 + j*numPoints + i] = (k == j) ? stabTurbEq : 0.0;

        if (computeJacobian == 1) {
            for (m = 0; m < nch; m++)
                for (k = 5; k < nch; k++)
                    for (j = 0; j < nch; j++)
                        tau_uh[m*sz3 + k*sz2 + j*numPoints + i] = (k == j) ? DstabTurbEq[m] : 0.0;
        }
    }
    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch

   // Make zero the drivatives for constant Roe:
   for (i=0; i<numPoints*nch*nch*nch; i++)
       tau_uh[i] = 0.0;

    delete[] DstabTurbEq;
}

void getRoeStabilizationTensor(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    if (nd == 2)
        getRoeStabilizationTensor2d(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (nd == 3)
        getRoeStabilizationTensor3d(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else
        error("Number of dimensions not implemented.\n");
    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch
}

void getProposedStabilizationTensor2d(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    // This functions works for compressible Euler, NS and RANS (any turbulence model).
    // If tauValue <= 0, the maximum eigenvalue is used to stabilize the turbulence model equation(s). Otherwise, tauValue is used.

    double tauValue;
    double nx, ny;
    double r, ru, rv, r1, uv, vv, un;
    double Vgx, Vgy, Vgn;

    double eps = 1.0e-8;        // eps > 0 ensures that tau is non-zero on at least one side of all faces even with finite-precision arithmetic

    tauValue = param[5];

    int sz2 = nch * numPoints;
    int i, j, k, kk;

    Int ALEflag = app.ALEflag;

    for (i=0; i<numPoints; i++) {
        r   = uh[0*numPoints+i]; // 1.0000;
        ru  = uh[1*numPoints+i]; // 0.9994;
        rv  = uh[2*numPoints+i]; // 0.0349;

        nx = nl[0 * numPoints + i];
        ny = nl[1 * numPoints + i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
        }

        r1 = 1.0 / r;
        uv = ru * r1;
        vv = rv * r1;
        un = uv * nx + vv * ny;

        Vgn = Vgx * nx + Vgy * ny;

        if (un - Vgn >= - eps) {
            for (k = 0; k < nch; k++)
                for (j = 0; j < nch; j++)
                    tau[i + j * numPoints + k * sz2] = (j == k) ? tauValue : 0.0;

            if (computeJacobian == 1) {
                for (kk = 0; kk < nch; kk++)
                    for (k = 0; k < nch; k++)
                        for (j = 0; j < nch; j++)
                            tau_uh[i + j * numPoints + k * numPoints * nch + kk * numPoints * nch * nch] = 0.0;
            }
        }
        else {
            for (k = 0; k < nch; k++)
                for (j = 0; j < nch; j++)
                    tau[i + j * numPoints + k * sz2] = (j == k) ? 0.0 : 0.0;

            if (computeJacobian == 1) {
                for (kk = 0; kk < nch; kk++)
                    for (k = 0; k < nch; k++)
                        for (j = 0; j < nch; j++)
                            tau_uh[i + j * numPoints + k * numPoints * nch + kk * numPoints * nch * nch] = 0.0;
            }
        }
    }
}


void getProposedStabilizationTensor3d(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    // This functions works for compressible Euler, NS and RANS (any turbulence model).
    // If tauValue <= 0, the maximum eigenvalue is used to stabilize the turbulence model equation(s). Otherwise, tauValue is used.

    double tauValue;
    double nx, ny, nz;
    double r, ru, rv, rw, r1, uv, vv, wv, un;
    double Vgx, Vgy, Vgz, Vgn;

    double eps = 1.0e-6;        // eps > 0 ensures that tau is non-zero on at least one side of all faces even with finite-precision arithmetic

    tauValue = param[5];

    int sz2 = nch * numPoints;
    int i, j, k, kk;

    Int ALEflag = app.ALEflag;

    for (i=0; i<numPoints; i++) {
        r   = uh[0*numPoints+i];
        ru  = uh[1*numPoints+i];
        rv  = uh[2*numPoints+i];
        rw  = uh[3*numPoints+i];

        nx = nl[0 * numPoints + i];
        ny = nl[1 * numPoints + i];
        nz = nl[2 * numPoints + i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
            Vgz = pg[(3 * nd + 1) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
            Vgz = 0.0;
        }

        r1 = 1.0 / r;
        uv = ru * r1;
        vv = rv * r1;
        wv = rw * r1;
        un = uv * nx + vv * ny + wv * nz;

        Vgn = Vgx * nx + Vgy * ny + Vgz * nz;

        if (un - Vgn >= - eps) {
            for (k = 0; k < nch; k++)
                for (j = 0; j < nch; j++)
                    tau[i + j * numPoints + k * sz2] = (j == k) ? tauValue : 0.0;

            if (computeJacobian == 1) {
                for (kk = 0; kk < nch; kk++)
                    for (k = 0; k < nch; k++)
                        for (j = 0; j < nch; j++)
                            tau_uh[i + j * numPoints + k * numPoints * nch + kk * numPoints * nch * nch] = 0.0;
            }
        }
        else {
            for (k = 0; k < nch; k++)
                for (j = 0; j < nch; j++)
                    tau[i + j * numPoints + k * sz2] = 0.0;

            if (computeJacobian == 1) {
                for (kk = 0; kk < nch; kk++)
                    for (k = 0; k < nch; k++)
                        for (j = 0; j < nch; j++)
                            tau_uh[i + j * numPoints + k * numPoints * nch + kk * numPoints * nch * nch] = 0.0;
            }
        }
    }
}


void getProposedStabilizationTensor(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    if (nd == 2)
        getProposedStabilizationTensor2d(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (nd == 3)
        getProposedStabilizationTensor3d(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else
        error("Number of dimensions not implemented.\n");

    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch
}


void getUpwindedStabilizationTensor2d(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    // This functions works for compressible Euler, NS and RANS (any turbulence model).
    // If tauValue <= 0, the maximum eigenvalue is used to stabilize the turbulence model equation(s). Otherwise, tauValue is used.

    double gam, gam1, epslm, tauValue, signun, ic, Dc2Reg_Dc2;
    double nx, ny;
    double r, ru, rv, rE;
    double r1, r1m1, r1m2, r1m3, r1m4;
    double uv, uvm1, uvm2, uvm3, uvm4;
    double vv, vvm1, vvm2, vvm3, vvm4;
    double Vgx, Vgy, Vgn;
    double E, Em1, Em2, Em3, Em4;
    double af, afm1, afm2, afm3, afm4;
    double p, pm1, pm2, pm3, pm4;
    double h, hm1, hm2, hm3, hm4;
    double s1, s1m1, s1m2, s1m3, s1m4;
    double s2, s2m1, s2m2, s2m3, s2m4;
    double c2, c2m1, c2m2, c2m3, c2m4;
    double c, cm1, cm2, cm3, cm4;
    double un, unm1, unm2, unm3, unm4;
    double cc1, cc1m1, cc1m2, cc1m3, cc1m4;
    double cc2, cc2m1, cc2m2, cc2m3, cc2m4;
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4;
    double rlam1, rlam1m1, rlam1m2, rlam1m3, rlam1m4;
    double rlam2, rlam2m1, rlam2m2, rlam2m3, rlam2m4;
    double rlam3, rlam3m1, rlam3m2, rlam3m3, rlam3m4;
    double lamMax, lamMaxm1, lamMaxm2, lamMaxm3, lamMaxm4;
    
    double x, Drlam1Dx, Drlam2Dx, Drlam3Dx;
    double eps = 1.0e-8;        // eps > 0 ensures that tau is non-zero on at least one side of all faces even with finite-precision arithmetic
    double pi = 3.141592653589793;
    double b = 100.0;
    double minTau = 0.1;
    
    gam   = param[0];
    gam1 = gam - 1.0;
    epslm = param[1];
    tauValue = param[5];
    
    double stabTurbEq;
    double *DstabTurbEq  = new double [nch];

    int nch2 = nch * nch, sz2 = nch * numPoints, sz3 = sz2 * nch;
    int i, j, k, m, index_tmp;
    
    Int ALEflag = app.ALEflag;
    
    for (i=0; i<numPoints; i++) {
       r   = uh[0*numPoints+i]; // 1.0000;
       ru  = uh[1*numPoints+i]; // 0.9994;
       rv  = uh[2*numPoints+i]; // 0.0349;
       rE  = uh[3*numPoints+i]; // 45.1429; // 1.7858e+04

        nx  = nl[0*numPoints+i];
        ny  = nl[1*numPoints+i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
        }

        r1   = 1.0/r;
        uv   = ru*r1;
        vv   = rv*r1;
        E    = rE*r1;
        af   = 0.5*(uv*uv+vv*vv);
        p    = gam1*(rE -r*af);
        h    = E   + p*r1;
        c2   = gam*p*r1;
        Dc2Reg_Dc2 = atan(b*c2)/pi + 1.0/2.0 + (b*c2) / ((1+b*c2*b*c2) * pi);
        c2 = max ( c2*(atan(b*c2)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 , 1.0e-6);
        c    = sqrt(c2);
        un   = uv*nx   + vv*ny;

        r1m1 = -1.0/(r*r);
        r1m2 = 0.0;
        r1m3 = 0.0;
        r1m4 = 0.0;

        if (computeJacobian == 1) {
            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;
            uvm3 = ru * r1m3;
            uvm4 = ru * r1m4;

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;
            vvm3 = r1 + rv * r1m3;
            vvm4 = rv * r1m4;

            Em1 = rE * r1m1;
            Em2 = rE * r1m2;
            Em3 = rE * r1m3;
            Em4 = r1 + rE * r1m4;

            afm1 = uv * uvm1 + vv * vvm1;
            afm2 = uv * uvm2 + vv * vvm2;
            afm3 = uv * uvm3 + vv * vvm3;
            afm4 = uv * uvm4 + vv * vvm4;

            pm1 = gam1 * (-af - r * afm1);
            pm2 = gam1 * (-r * afm2);
            pm3 = gam1 * (-r * afm3);
            pm4 = gam1 * (1.0 - r * afm4);

            hm1 = Em1 + pm1 * r1 + p * r1m1;
            hm2 = Em2 + pm2 * r1 + p * r1m2;
            hm3 = Em3 + pm3 * r1 + p * r1m3;
            hm4 = Em4 + pm4 * r1 + p * r1m4;

            c2m1 = gam * (pm1 * r1 + p * r1m1);
            c2m2 = gam * (pm2 * r1 + p * r1m2);
            c2m3 = gam * (pm3 * r1 + p * r1m3);
            c2m4 = gam * (pm4 * r1 + p * r1m4);
            
            c2m1 *= Dc2Reg_Dc2;
            c2m2 *= Dc2Reg_Dc2;
            c2m3 *= Dc2Reg_Dc2;
            c2m4 *= Dc2Reg_Dc2;
            
            cm1 = 0.5 * c2m1 / c;
            cm2 = 0.5 * c2m2 / c;
            cm3 = 0.5 * c2m3 / c;
            cm4 = 0.5 * c2m4 / c;

            unm1 = uvm1 * nx + vvm1 * ny;
            unm2 = uvm2 * nx + vvm2 * ny;
            unm3 = uvm3 * nx + vvm3 * ny;
            unm4 = uvm4 * nx + vvm4 * ny;
        }

        Vgn = Vgx * nx + Vgy * ny;

        x = un - Vgn + c - minTau;
        rlam1 = x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 + minTau;
        Drlam1Dx = atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi);
        signun = (un - Vgn + c < 0) ? -1.0 : 1.0;
        if (computeJacobian == 1) {
            rlam1m1 = Drlam1Dx * (unm1 + cm1);
            rlam1m2 = Drlam1Dx * (unm2 + cm2);
            rlam1m3 = Drlam1Dx * (unm3 + cm3);
            rlam1m4 = Drlam1Dx * (unm4 + cm4);
        }

        if (epslm > 0 && (un - Vgn + c > 0)) {
            rlam = 0.5 * (rlam1 * rlam1 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam1 * rlam1m1 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam1 * rlam1m2 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam1 * rlam1m3 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam1 * rlam1m4 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm4);
            }
        
            ic = rlam1 < epslm * c;
            rlam1 = ic * rlam + (1 - ic) * rlam1;
            
            if (computeJacobian == 1) {
                rlam1m1 = ic * rlamm1 + (1 - ic) * rlam1m1;
                rlam1m2 = ic * rlamm2 + (1 - ic) * rlam1m2;
                rlam1m3 = ic * rlamm3 + (1 - ic) * rlam1m3;
                rlam1m4 = ic * rlamm4 + (1 - ic) * rlam1m4;
            }
        }

        x = un - Vgn - c - minTau;
        rlam2 = x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 + minTau;
        Drlam2Dx = atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi);
        signun = (un - Vgn - c < 0) ? -1.0 : 1.0;
        if (computeJacobian == 1) {
            rlam2m1 = Drlam2Dx * (unm1 - cm1);
            rlam2m2 = Drlam2Dx * (unm2 - cm2);
            rlam2m3 = Drlam2Dx * (unm3 - cm3);
            rlam2m4 = Drlam2Dx * (unm4 - cm4);
        }

        if (epslm > 0 && (un - Vgn - c > 0)) {
            rlam = 0.5 * (rlam2 * rlam2 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam2 * rlam2m1 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam2 * rlam2m2 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam2 * rlam2m3 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam2 * rlam2m4 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm4);
            }

            ic = rlam2 < epslm * c;
            rlam2 = ic * rlam + (1 - ic) * rlam2;

            if (computeJacobian == 1) {
                rlam2m1 = ic * rlamm1 + (1 - ic) * rlam2m1;
                rlam2m2 = ic * rlamm2 + (1 - ic) * rlam2m2;
                rlam2m3 = ic * rlamm3 + (1 - ic) * rlam2m3;
                rlam2m4 = ic * rlamm4 + (1 - ic) * rlam2m4;
            }
        }

        x = un - Vgn - minTau;
        rlam3 = x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 + minTau;
        Drlam3Dx = atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi);
        signun = (un - Vgn < 0) ? -1.0 : 1.0;
        if (computeJacobian == 1) {
            rlam3m1 = Drlam3Dx * unm1;
            rlam3m2 = Drlam3Dx * unm2;
            rlam3m3 = Drlam3Dx * unm3;
            rlam3m4 = Drlam3Dx * unm4;
        }

        if (epslm > 0 && (un - Vgn > 0)) {
            rlam = 0.5 * (rlam3 * rlam3 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam3 * rlam3m1 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam3 * rlam3m2 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam3 * rlam3m3 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam3 * rlam3m4 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm4);
            }

            ic = rlam3 < epslm * c;
            rlam3 = ic * rlam + (1 - ic) * rlam3;

            if (computeJacobian == 1) {
                rlam3m1 = ic * rlamm1 + (1 - ic) * rlam3m1;
                rlam3m2 = ic * rlamm2 + (1 - ic) * rlam3m2;
                rlam3m3 = ic * rlamm3 + (1 - ic) * rlam3m3;
                rlam3m4 = ic * rlamm4 + (1 - ic) * rlam3m4;
            }
        }

        if (un - Vgn >= 0) {
            lamMax = (un-Vgn) + c;
            if (computeJacobian == 1) {
                lamMaxm1 = unm1 + cm1;
                lamMaxm2 = unm2 + cm2;
                lamMaxm3 = unm3 + cm3;
                lamMaxm4 = unm4 + cm4;
            }
        }
        else {
            lamMax = - (un-Vgn) + c;
            if (computeJacobian == 1) {
                lamMaxm1 = -unm1 + cm1;
                lamMaxm2 = -unm2 + cm2;
                lamMaxm3 = -unm3 + cm3;
                lamMaxm4 = -unm4 + cm4;
            }
        }

        if (tauValue > 0) {
            stabTurbEq = tauValue;
            if (computeJacobian == 1) {
                for (j = 0; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }
        else {
            stabTurbEq = lamMax;
            if (computeJacobian == 1) {
                DstabTurbEq[0] = lamMaxm1;
                DstabTurbEq[1] = lamMaxm2;
                DstabTurbEq[2] = lamMaxm3;
                DstabTurbEq[3] = lamMaxm4;
                for (j = 4; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }

        s1      = 0.5*(rlam1   + rlam2);
        s2      = 0.5*(rlam1   - rlam2);

        if (computeJacobian == 1) {
            s1m1 = 0.5 * (rlam1m1 + rlam2m1);
            s1m2 = 0.5 * (rlam1m2 + rlam2m2);
            s1m3 = 0.5 * (rlam1m3 + rlam2m3);
            s1m4 = 0.5 * (rlam1m4 + rlam2m4);

            s2m1 = 0.5 * (rlam1m1 - rlam2m1);
            s2m2 = 0.5 * (rlam1m2 - rlam2m2);
            s2m3 = 0.5 * (rlam1m3 - rlam2m3);
            s2m4 = 0.5 * (rlam1m4 - rlam2m4);
        }

        cc1   = gam1*(s1-rlam3)*af/c2-(s2*un/c);
        cc2   = gam1*s2*af/c-(s1-rlam3)*un;

        index_tmp = 0*nch*numPoints+i;
        tau[index_tmp + 0*numPoints]  = rlam3+cc1;
        tau[index_tmp + 1*numPoints]  = cc1*uv+cc2*nx;
        tau[index_tmp + 2*numPoints]  = cc1*vv+cc2*ny;
        tau[index_tmp + 3*numPoints]  = cc1*h+cc2*un;
        for (j=4; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = gam1*((s1m1-rlam3m1)*af/c2 + (s1-rlam3)*afm1/c2 - (s1-rlam3)*af*c2m1/(c2*c2)) - s2m1*un/c - s2*unm1/c + s2*un*cm1/(c*c);
            cc1m2 = gam1*((s1m2-rlam3m2)*af/c2 + (s1-rlam3)*afm2/c2 - (s1-rlam3)*af*c2m2/(c2*c2)) - s2m2*un/c - s2*unm2/c + s2*un*cm2/(c*c);
            cc1m3 = gam1*((s1m3-rlam3m3)*af/c2 + (s1-rlam3)*afm3/c2 - (s1-rlam3)*af*c2m3/(c2*c2)) - s2m3*un/c - s2*unm3/c + s2*un*cm3/(c*c);
            cc1m4 = gam1*((s1m4-rlam3m4)*af/c2 + (s1-rlam3)*afm4/c2 - (s1-rlam3)*af*c2m4/(c2*c2)) - s2m4*un/c - s2*unm4/c + s2*un*cm4/(c*c);
            
            cc2m1 = gam1*(s2m1*af/c + s2*afm1/c - s2*af*cm1/(c*c)) - (s1m1-rlam3m1)*un - (s1-rlam3)*unm1;
            cc2m2 = gam1*(s2m2*af/c + s2*afm2/c - s2*af*cm2/(c*c)) - (s1m2-rlam3m2)*un - (s1-rlam3)*unm2;
            cc2m3 = gam1*(s2m3*af/c + s2*afm3/c - s2*af*cm3/(c*c)) - (s1m3-rlam3m3)*un - (s1-rlam3)*unm3;
            cc2m4 = gam1*(s2m4*af/c + s2*afm4/c - s2*af*cm4/(c*c)) - (s1m4-rlam3m4)*un - (s1-rlam3)*unm4;

            index_tmp = (0*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m1 + cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m2 + cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m3 + cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m4 + cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            
            index_tmp = 0*sz2 + i;
            for (k=4; k<nch; k++)
                for (j=0; j<4; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=4; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }
        
        cc1   = -gam1*(s1-rlam3)*uv/c2+(s2*nx/c);
        cc2   = -gam1*s2*uv/c + (s1-rlam3)*nx;
        
        index_tmp = 1*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = rlam3+cc1*uv+cc2*nx;
        tau[index_tmp + 2*numPoints] = cc1*vv+cc2*ny;
        tau[index_tmp + 3*numPoints] = cc1*h+cc2*un;
        for (j=4; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = -gam1*((s1m1-rlam3m1)*uv/c2 + (s1-rlam3)*uvm1/c2 - (s1-rlam3)*uv*c2m1/(c2*c2)) + s2m1*nx/c - s2*nx*cm1/(c*c);
            cc1m2 = -gam1*((s1m2-rlam3m2)*uv/c2 + (s1-rlam3)*uvm2/c2 - (s1-rlam3)*uv*c2m2/(c2*c2)) + s2m2*nx/c - s2*nx*cm2/(c*c);
            cc1m3 = -gam1*((s1m3-rlam3m3)*uv/c2 + (s1-rlam3)*uvm3/c2 - (s1-rlam3)*uv*c2m3/(c2*c2)) + s2m3*nx/c - s2*nx*cm3/(c*c);
            cc1m4 = -gam1*((s1m4-rlam3m4)*uv/c2 + (s1-rlam3)*uvm4/c2 - (s1-rlam3)*uv*c2m4/(c2*c2)) + s2m4*nx/c - s2*nx*cm4/(c*c);

            cc2m1 = -gam1*(s2m1*uv/c + s2*uvm1/c - s2*uv*cm1/(c*c)) + (s1m1-rlam3m1)*nx;
            cc2m2 = -gam1*(s2m2*uv/c + s2*uvm2/c - s2*uv*cm2/(c*c)) + (s1m2-rlam3m2)*nx;
            cc2m3 = -gam1*(s2m3*uv/c + s2*uvm3/c - s2*uv*cm3/(c*c)) + (s1m3-rlam3m3)*nx;
            cc2m4 = -gam1*(s2m4*uv/c + s2*uvm4/c - s2*uv*cm4/(c*c)) + (s1m4-rlam3m4)*nx;
            
            index_tmp = (0*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = rlam3m1 + cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = rlam3m2 + cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = rlam3m3 + cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = rlam3m4 + cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;

            index_tmp = 1*sz2 + i;
            for (k=4; k<nch; k++)
                for (j=0; j<4; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=4; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }
        
        cc1   = -gam1*(s1-rlam3)*vv/c2+(s2*ny/c);
        cc2   = -gam1*s2*vv/c+(s1-rlam3)*ny;

        index_tmp = 2*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = cc1*uv+cc2*nx;
        tau[index_tmp + 2*numPoints] = rlam3+cc1*vv+cc2*ny;
        tau[index_tmp + 3*numPoints] = cc1*h+cc2*un;
        for (j=4; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = -gam1*((s1m1-rlam3m1)*vv/c2 + (s1-rlam3)*vvm1/c2 - (s1-rlam3)*vv*c2m1/(c2*c2)) + s2m1*ny/c - s2*ny*cm1/(c*c);
            cc1m2 = -gam1*((s1m2-rlam3m2)*vv/c2 + (s1-rlam3)*vvm2/c2 - (s1-rlam3)*vv*c2m2/(c2*c2)) + s2m2*ny/c - s2*ny*cm2/(c*c);
            cc1m3 = -gam1*((s1m3-rlam3m3)*vv/c2 + (s1-rlam3)*vvm3/c2 - (s1-rlam3)*vv*c2m3/(c2*c2)) + s2m3*ny/c - s2*ny*cm3/(c*c);
            cc1m4 = -gam1*((s1m4-rlam3m4)*vv/c2 + (s1-rlam3)*vvm4/c2 - (s1-rlam3)*vv*c2m4/(c2*c2)) + s2m4*ny/c - s2*ny*cm4/(c*c);

            cc2m1 = -gam1*(s2m1*vv/c + s2*vvm1/c - s2*vv*cm1/(c*c)) + (s1m1-rlam3m1)*ny;
            cc2m2 = -gam1*(s2m2*vv/c + s2*vvm2/c - s2*vv*cm2/(c*c)) + (s1m2-rlam3m2)*ny;
            cc2m3 = -gam1*(s2m3*vv/c + s2*vvm3/c - s2*vv*cm3/(c*c)) + (s1m3-rlam3m3)*ny;
            cc2m4 = -gam1*(s2m4*vv/c + s2*vvm4/c - s2*vv*cm4/(c*c)) + (s1m4-rlam3m4)*ny;

            index_tmp = (0*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m1 + cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m2 + cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m3 + cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m4 + cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;

            index_tmp = 2*sz2 + i;
            for (k=4; k<nch; k++)
                for (j=0; j<4; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=4; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1   = gam1*(s1-rlam3)/c2;
        cc2   = gam1*s2/c;

        index_tmp = 3*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = cc1*uv+cc2*nx;
        tau[index_tmp + 2*numPoints] = cc1*vv+cc2*ny;
        tau[index_tmp + 3*numPoints] = rlam3+cc1*h+cc2*un;
        for (j=4; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = gam1*((s1m1-rlam3m1)/c2 - (s1-rlam3)*c2m1/(c2*c2));
            cc1m2 = gam1*((s1m2-rlam3m2)/c2 - (s1-rlam3)*c2m2/(c2*c2));
            cc1m3 = gam1*((s1m3-rlam3m3)/c2 - (s1-rlam3)*c2m3/(c2*c2));
            cc1m4 = gam1*((s1m4-rlam3m4)/c2 - (s1-rlam3)*c2m4/(c2*c2));

            cc2m1 = gam1*(s2m1/c - s2*cm1/(c*c));
            cc2m2 = gam1*(s2m2/c - s2*cm2/(c*c));
            cc2m3 = gam1*(s2m3/c - s2*cm3/(c*c));
            cc2m4 = gam1*(s2m4/c - s2*cm4/(c*c));

            index_tmp = (0*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m1 + cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m2 + cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m3 + cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m4 + cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;

            index_tmp = 3*sz2 + i;
            for (k=4; k<nch; k++)
                for (j=0; j<4; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=4; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        for (k=4; k<nch; k++)
            for (j=0; j<nch; j++)
                tau[k*sz2 + j*numPoints + i] = (k == j) ? stabTurbEq : 0.0;

        if (computeJacobian == 1) {
            for (m = 0; m < nch; m++)
                for (k = 4; k < nch; k++)
                    for (j = 0; j < nch; j++)
                        tau_uh[m*sz3 + k*sz2 + j*numPoints + i] = (k == j) ? DstabTurbEq[m] : 0.0;
        }
    }
    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch
    
    for (i=0; i<numPoints*nch*nch*nch; i++)
        tau_uh[i] = 0.0;

    delete[] DstabTurbEq;
}

void getUpwindedStabilizationTensor3d(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    // This functions works for compressible Euler, NS and RANS (any turbulence model).
    // If tauValue <= 0, the maximum eigenvalue is used to stabilize the turbulence model equation(s). Otherwise, tauValue is used.

//    getan3d(tau, tau_uh, uh, nl, param, 0, numPoints);

    double gam, gam1, epslm, tauValue, signun, ic, Dc2Reg_Dc2;
    double nx, ny, nz;
    double r, ru, rv, rw, rE;
    double r1, r1m1, r1m2, r1m3, r1m4, r1m5;
    double uv, uvm1, uvm2, uvm3, uvm4, uvm5;
    double vv, vvm1, vvm2, vvm3, vvm4, vvm5;
    double wv, wvm1, wvm2, wvm3, wvm4, wvm5;
    double Vgx, Vgy, Vgz, Vgn;
    double E, Em1, Em2, Em3, Em4, Em5;
    double af, afm1, afm2, afm3, afm4, afm5;
    double p, pm1, pm2, pm3, pm4, pm5;
    double h, hm1, hm2, hm3, hm4, hm5;
    double s1, s1m1, s1m2, s1m3, s1m4, s1m5;
    double s2, s2m1, s2m2, s2m3, s2m4, s2m5;
    double c2, c2m1, c2m2, c2m3, c2m4, c2m5;
    double c, cm1, cm2, cm3, cm4, cm5;
    double un, unm1, unm2, unm3, unm4, unm5;
    double cc1, cc1m1, cc1m2, cc1m3, cc1m4, cc1m5;
    double cc2, cc2m1, cc2m2, cc2m3, cc2m4, cc2m5;
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4, rlamm5;
    double rlam1, rlam1m1, rlam1m2, rlam1m3, rlam1m4, rlam1m5;
    double rlam2, rlam2m1, rlam2m2, rlam2m3, rlam2m4, rlam2m5;
    double rlam3, rlam3m1, rlam3m2, rlam3m3, rlam3m4, rlam3m5;
    double lamMax, lamMaxm1, lamMaxm2, lamMaxm3, lamMaxm4, lamMaxm5;

    double x, Drlam1Dx, Drlam2Dx, Drlam3Dx;
    double eps = 1.0e-8;        // eps > 0 ensures that tau is non-zero on at least one side of all faces even with finite-precision arithmetic
    double pi = 3.141592653589793;
    double b = 1.0e4;
    double minTau = 0.1;
    
    gam = param[0];
    gam1 = gam - 1.0;
    epslm = param[1];
    tauValue = param[5];

    double stabTurbEq;
    double *DstabTurbEq  = new double [nch];

    int nch2 = nch * nch, sz2 = nch * numPoints, sz3 = sz2 * nch;
    int i, j, k, m, index_tmp;

    Int ALEflag = app.ALEflag;

    for (i = 0; i < numPoints; i++) {
        nx = nl[0 * numPoints + i];
        ny = nl[1 * numPoints + i];
        nz = nl[2 * numPoints + i];

        r = uh[0 * numPoints + i];
        ru = uh[1 * numPoints + i];
        rv = uh[2 * numPoints + i];
        rw = uh[3 * numPoints + i];
        rE = uh[4 * numPoints + i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
            Vgz = pg[(2 * nd + 2) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
            Vgz = 0.0;
        }

        r1 = 1.0 / r;
        uv = ru * r1;
        vv = rv * r1;
        wv = rw * r1;
        E = rE * r1;
        af = 0.5 * (uv * uv + vv * vv + wv * wv);
        p = gam1 * (rE - r * af);
        h = E + p * r1;
        c2 = gam * p * r1;
        Dc2Reg_Dc2 = atan(b*c2)/pi + 1.0/2.0 + (b*c2) / ((1+b*c2*b*c2) * pi);
        c2 = max ( c2*(atan(b*c2)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 , 1.0e-6);
        c = sqrt(c2);
        un = uv * nx + vv * ny + wv * nz;

        if (computeJacobian == 1) {
            r1m1 = -1.0 / (r * r);
            r1m2 = 0.0;
            r1m3 = 0.0;
            r1m4 = 0.0;
            r1m5 = 0.0;

            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;
            uvm3 = ru * r1m3;
            uvm4 = ru * r1m4;
            uvm5 = ru * r1m5;

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;
            vvm3 = r1 + rv * r1m3;
            vvm4 = rv * r1m4;
            vvm5 = rv * r1m5;

            wvm1 = rw * r1m1;
            wvm2 = rw * r1m2;
            wvm3 = rw * r1m3;
            wvm4 = r1 + rw * r1m4;
            wvm5 = rw * r1m5;

            Em1 = rE * r1m1;
            Em2 = rE * r1m2;
            Em3 = rE * r1m3;
            Em4 = rE * r1m4;
            Em5 = r1 + rE * r1m5;

            afm1 = uv * uvm1 + vv * vvm1 + wv * wvm1;
            afm2 = uv * uvm2 + vv * vvm2 + wv * wvm2;
            afm3 = uv * uvm3 + vv * vvm3 + wv * wvm3;
            afm4 = uv * uvm4 + vv * vvm4 + wv * wvm4;
            afm5 = uv * uvm5 + vv * vvm5 + wv * wvm5;

            pm1 = gam1 * (-af - r * afm1);
            pm2 = gam1 * (-r * afm2);
            pm3 = gam1 * (-r * afm3);
            pm4 = gam1 * (-r * afm4);
            pm5 = gam1 * (1.0 - r * afm5);

            hm1 = Em1 + pm1 * r1 + p * r1m1;
            hm2 = Em2 + pm2 * r1 + p * r1m2;
            hm3 = Em3 + pm3 * r1 + p * r1m3;
            hm4 = Em4 + pm4 * r1 + p * r1m4;
            hm5 = Em5 + pm5 * r1 + p * r1m5;

            c2m1 = gam * (pm1 * r1 + p * r1m1);
            c2m2 = gam * (pm2 * r1 + p * r1m2);
            c2m3 = gam * (pm3 * r1 + p * r1m3);
            c2m4 = gam * (pm4 * r1 + p * r1m4);
            c2m5 = gam * (pm5 * r1 + p * r1m5);
            
            c2m1 *= Dc2Reg_Dc2;
            c2m2 *= Dc2Reg_Dc2;
            c2m3 *= Dc2Reg_Dc2;
            c2m4 *= Dc2Reg_Dc2;
            c2m5 *= Dc2Reg_Dc2;

            cm1 = 0.5 * c2m1 / c;
            cm2 = 0.5 * c2m2 / c;
            cm3 = 0.5 * c2m3 / c;
            cm4 = 0.5 * c2m4 / c;
            cm5 = 0.5 * c2m5 / c;

            unm1 = uvm1 * nx + vvm1 * ny + wvm1 * nz;
            unm2 = uvm2 * nx + vvm2 * ny + wvm2 * nz;
            unm3 = uvm3 * nx + vvm3 * ny + wvm3 * nz;
            unm4 = uvm4 * nx + vvm4 * ny + wvm4 * nz;
            unm5 = uvm5 * nx + vvm5 * ny + wvm5 * nz;
        }

        Vgn = Vgx * nx + Vgy * ny + Vgz * nz;
        
        x = un - Vgn + c - minTau;
        rlam1 = x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 + minTau;
        Drlam1Dx = atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi);
        signun = (un - Vgn + c < 0) ? -1.0 : 1.0;
        if (computeJacobian == 1) {
            rlam1m1 = Drlam1Dx * (unm1 + cm1);
            rlam1m2 = Drlam1Dx * (unm2 + cm2);
            rlam1m3 = Drlam1Dx * (unm3 + cm3);
            rlam1m4 = Drlam1Dx * (unm4 + cm4);
            rlam1m5 = Drlam1Dx * (unm5 + cm5);
        }

        if (epslm > 0 && (un - Vgn + c > 0)) {
            rlam = 0.5 * (rlam1 * rlam1 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam1 * rlam1m1 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam1 * rlam1m2 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam1 * rlam1m3 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam1 * rlam1m4 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm4);
                rlamm5 =
                        rlam1 * rlam1m5 / (epslm * c) + 0.5 * (1 - rlam1 * rlam1 / (epslm * c * epslm * c)) * (epslm * cm5);
            }

            ic = rlam1 < epslm * c;
            rlam1 = ic * rlam + (1 - ic) * rlam1;

            if (computeJacobian == 1) {
                rlam1m1 = ic * rlamm1 + (1 - ic) * rlam1m1;
                rlam1m2 = ic * rlamm2 + (1 - ic) * rlam1m2;
                rlam1m3 = ic * rlamm3 + (1 - ic) * rlam1m3;
                rlam1m4 = ic * rlamm4 + (1 - ic) * rlam1m4;
                rlam1m5 = ic * rlamm5 + (1 - ic) * rlam1m5;
            }
        }

        x = un - Vgn - c - minTau;
        rlam2 = x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 + minTau;
        Drlam2Dx = atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi);
        signun = (un - Vgn - c < 0) ? -1.0 : 1.0;
        if (computeJacobian == 1) {
            rlam2m1 = Drlam2Dx * (unm1 - cm1);
            rlam2m2 = Drlam2Dx * (unm2 - cm2);
            rlam2m3 = Drlam2Dx * (unm3 - cm3);
            rlam2m4 = Drlam2Dx * (unm4 - cm4);
            rlam2m5 = Drlam2Dx * (unm5 - cm5);
        }
        
        if (epslm > 0 && (un - Vgn - c > 0)) {
            rlam = 0.5 * (rlam2 * rlam2 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam2 * rlam2m1 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam2 * rlam2m2 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam2 * rlam2m3 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam2 * rlam2m4 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm4);
                rlamm5 =
                        rlam2 * rlam2m5 / (epslm * c) + 0.5 * (1 - rlam2 * rlam2 / (epslm * c * epslm * c)) * (epslm * cm5);
            }

            ic = rlam2 < epslm * c;
            rlam2 = ic * rlam + (1 - ic) * rlam2;

            if (computeJacobian == 1) {
                rlam2m1 = ic * rlamm1 + (1 - ic) * rlam2m1;
                rlam2m2 = ic * rlamm2 + (1 - ic) * rlam2m2;
                rlam2m3 = ic * rlamm3 + (1 - ic) * rlam2m3;
                rlam2m4 = ic * rlamm4 + (1 - ic) * rlam2m4;
                rlam2m5 = ic * rlamm5 + (1 - ic) * rlam2m5;
            }
        }

        x = un - Vgn - minTau;
        rlam3 = x*(atan(b*x)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 + minTau;
        Drlam3Dx = atan(b*x)/pi + 1.0/2.0 + (b*x) / ((1+b*x*b*x) * pi);
        signun = (un - Vgn < 0) ? -1.0 : 1.0;
        if (computeJacobian == 1) {
            rlam3m1 = Drlam3Dx * unm1;
            rlam3m2 = Drlam3Dx * unm2;
            rlam3m3 = Drlam3Dx * unm3;
            rlam3m4 = Drlam3Dx * unm4;
            rlam3m5 = Drlam3Dx * unm5;
        }
        
        if (epslm > 0 && (un - Vgn > 0)) {
            rlam = 0.5 * (rlam3 * rlam3 / (epslm * c) + epslm * c);

            if (computeJacobian == 1) {
                rlamm1 =
                        rlam3 * rlam3m1 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm1);
                rlamm2 =
                        rlam3 * rlam3m2 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm2);
                rlamm3 =
                        rlam3 * rlam3m3 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm3);
                rlamm4 =
                        rlam3 * rlam3m4 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm4);
                rlamm5 =
                        rlam3 * rlam3m5 / (epslm * c) + 0.5 * (1 - rlam3 * rlam3 / (epslm * c * epslm * c)) * (epslm * cm5);
            }

            ic = rlam3 < epslm * c;
            rlam3 = ic * rlam + (1 - ic) * rlam3;

            if (computeJacobian == 1) {
                rlam3m1 = ic * rlamm1 + (1 - ic) * rlam3m1;
                rlam3m2 = ic * rlamm2 + (1 - ic) * rlam3m2;
                rlam3m3 = ic * rlamm3 + (1 - ic) * rlam3m3;
                rlam3m4 = ic * rlamm4 + (1 - ic) * rlam3m4;
                rlam3m5 = ic * rlamm5 + (1 - ic) * rlam3m5;
            }
        }

        if (un - Vgn >= 0) {
            lamMax = (un-Vgn) + c;
            if (computeJacobian == 1) {
                lamMaxm1 = unm1 + cm1;
                lamMaxm2 = unm2 + cm2;
                lamMaxm3 = unm3 + cm3;
                lamMaxm4 = unm4 + cm4;
                lamMaxm5 = unm5 + cm5;
            }
        }
        else {
            lamMax = - (un-Vgn) + c;
            if (computeJacobian == 1) {
                lamMaxm1 = -unm1 + cm1;
                lamMaxm2 = -unm2 + cm2;
                lamMaxm3 = -unm3 + cm3;
                lamMaxm4 = -unm4 + cm4;
                lamMaxm5 = -unm5 + cm5;
            }
        }

        if (tauValue > 0) {
            stabTurbEq = tauValue;
            if (computeJacobian == 1) {
                for (j = 0; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }
        else {
            stabTurbEq = lamMax;
            if (computeJacobian == 1) {
                DstabTurbEq[0] = lamMaxm1;
                DstabTurbEq[1] = lamMaxm2;
                DstabTurbEq[2] = lamMaxm3;
                DstabTurbEq[3] = lamMaxm4;
                DstabTurbEq[4] = lamMaxm5;
                for (j = 5; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }

        s1 = 0.5 * (rlam1 + rlam2);
        s2 = 0.5 * (rlam1 - rlam2);

        if (computeJacobian == 1) {
            s1m1 = 0.5 * (rlam1m1 + rlam2m1);
            s1m2 = 0.5 * (rlam1m2 + rlam2m2);
            s1m3 = 0.5 * (rlam1m3 + rlam2m3);
            s1m4 = 0.5 * (rlam1m4 + rlam2m4);
            s1m5 = 0.5 * (rlam1m5 + rlam2m5);

            s2m1 = 0.5 * (rlam1m1 - rlam2m1);
            s2m2 = 0.5 * (rlam1m2 - rlam2m2);
            s2m3 = 0.5 * (rlam1m3 - rlam2m3);
            s2m4 = 0.5 * (rlam1m4 - rlam2m4);
            s2m5 = 0.5 * (rlam1m5 - rlam2m5);
        }
        
        cc1 = gam1 * (s1 - rlam3) * af / c2 - (s2 * un / c);
        cc2 = gam1 * s2 * af / c - (s1 - rlam3) * un;

        index_tmp = 0*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = rlam3 + cc1;
        tau[index_tmp + 1*numPoints] = cc1 * uv + cc2 * nx;
        tau[index_tmp + 2*numPoints] = cc1 * vv + cc2 * ny;
        tau[index_tmp + 3*numPoints] = cc1 * wv + cc2 * nz;
        tau[index_tmp + 4*numPoints] = cc1 * h + cc2 * un;
        for (j=5; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = gam1 * ((s1m1 - rlam3m1) * af / c2 + (s1 - rlam3) * afm1 / c2 - (s1 - rlam3) * af * c2m1 / (c2 * c2)) -
                    s2m1 * un / c - s2 * unm1 / c + s2 * un * cm1 / (c * c);
            cc1m2 = gam1 * ((s1m2 - rlam3m2) * af / c2 + (s1 - rlam3) * afm2 / c2 - (s1 - rlam3) * af * c2m2 / (c2 * c2)) -
                    s2m2 * un / c - s2 * unm2 / c + s2 * un * cm2 / (c * c);
            cc1m3 = gam1 * ((s1m3 - rlam3m3) * af / c2 + (s1 - rlam3) * afm3 / c2 - (s1 - rlam3) * af * c2m3 / (c2 * c2)) -
                    s2m3 * un / c - s2 * unm3 / c + s2 * un * cm3 / (c * c);
            cc1m4 = gam1 * ((s1m4 - rlam3m4) * af / c2 + (s1 - rlam3) * afm4 / c2 - (s1 - rlam3) * af * c2m4 / (c2 * c2)) -
                    s2m4 * un / c - s2 * unm4 / c + s2 * un * cm4 / (c * c);
            cc1m5 = gam1 * ((s1m5 - rlam3m5) * af / c2 + (s1 - rlam3) * afm5 / c2 - (s1 - rlam3) * af * c2m5 / (c2 * c2)) -
                    s2m5 * un / c - s2 * unm5 / c + s2 * un * cm5 / (c * c);

            cc2m1 = gam1 * (s2m1 * af / c + s2 * afm1 / c - s2 * af * cm1 / (c * c)) - (s1m1 - rlam3m1) * un -
                    (s1 - rlam3) * unm1;
            cc2m2 = gam1 * (s2m2 * af / c + s2 * afm2 / c - s2 * af * cm2 / (c * c)) - (s1m2 - rlam3m2) * un -
                    (s1 - rlam3) * unm2;
            cc2m3 = gam1 * (s2m3 * af / c + s2 * afm3 / c - s2 * af * cm3 / (c * c)) - (s1m3 - rlam3m3) * un -
                    (s1 - rlam3) * unm3;
            cc2m4 = gam1 * (s2m4 * af / c + s2 * afm4 / c - s2 * af * cm4 / (c * c)) - (s1m4 - rlam3m4) * un -
                    (s1 - rlam3) * unm4;
            cc2m5 = gam1 * (s2m5 * af / c + s2 * afm5 / c - s2 * af * cm5 / (c * c)) - (s1m5 - rlam3m5) * un -
                    (s1 - rlam3) * unm5;

            index_tmp = (0*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m1 + cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m2 + cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m3 + cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m4 + cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            index_tmp += sz3; // (4*nch2 + 0*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = rlam3m5 + cc1m5;
            tau_uh[index_tmp + 1*numPoints] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;

            index_tmp = 0*sz2 + i;
            for (k=5; k<nch; k++)
                for (j=0; j<5; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=5; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1 = -gam1 * (s1 - rlam3) * uv / c2 + (s2 * nx / c);
        cc2 = -gam1 * s2 * uv / c + (s1 - rlam3) * nx;

        index_tmp = 1*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = rlam3 + cc1 * uv + cc2 * nx;
        tau[index_tmp + 2*numPoints] = cc1 * vv + cc2 * ny;
        tau[index_tmp + 3*numPoints] = cc1 * wv + cc2 * nz;
        tau[index_tmp + 4*numPoints] = cc1 * h + cc2 * un;
        for (j=5; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = -gam1 * ((s1m1 - rlam3m1) * uv / c2 + (s1 - rlam3) * uvm1 / c2 - (s1 - rlam3) * uv * c2m1 / (c2 * c2)) +
                    s2m1 * nx / c - s2 * nx * cm1 / (c * c);
            cc1m2 = -gam1 * ((s1m2 - rlam3m2) * uv / c2 + (s1 - rlam3) * uvm2 / c2 - (s1 - rlam3) * uv * c2m2 / (c2 * c2)) +
                    s2m2 * nx / c - s2 * nx * cm2 / (c * c);
            cc1m3 = -gam1 * ((s1m3 - rlam3m3) * uv / c2 + (s1 - rlam3) * uvm3 / c2 - (s1 - rlam3) * uv * c2m3 / (c2 * c2)) +
                    s2m3 * nx / c - s2 * nx * cm3 / (c * c);
            cc1m4 = -gam1 * ((s1m4 - rlam3m4) * uv / c2 + (s1 - rlam3) * uvm4 / c2 - (s1 - rlam3) * uv * c2m4 / (c2 * c2)) +
                    s2m4 * nx / c - s2 * nx * cm4 / (c * c);
            cc1m5 = -gam1 * ((s1m5 - rlam3m5) * uv / c2 + (s1 - rlam3) * uvm5 / c2 - (s1 - rlam3) * uv * c2m5 / (c2 * c2)) +
                    s2m5 * nx / c - s2 * nx * cm5 / (c * c);

            cc2m1 = -gam1 * (s2m1 * uv / c + s2 * uvm1 / c - s2 * uv * cm1 / (c * c)) + (s1m1 - rlam3m1) * nx;
            cc2m2 = -gam1 * (s2m2 * uv / c + s2 * uvm2 / c - s2 * uv * cm2 / (c * c)) + (s1m2 - rlam3m2) * nx;
            cc2m3 = -gam1 * (s2m3 * uv / c + s2 * uvm3 / c - s2 * uv * cm3 / (c * c)) + (s1m3 - rlam3m3) * nx;
            cc2m4 = -gam1 * (s2m4 * uv / c + s2 * uvm4 / c - s2 * uv * cm4 / (c * c)) + (s1m4 - rlam3m4) * nx;
            cc2m5 = -gam1 * (s2m5 * uv / c + s2 * uvm5 / c - s2 * uv * cm5 / (c * c)) + (s1m5 - rlam3m5) * nx;

            index_tmp = (0*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = rlam3m1 + cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = rlam3m2 + cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = rlam3m3 + cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = rlam3m4 + cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            index_tmp += sz3; // (4*nch2 + 1*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m5;
            tau_uh[index_tmp + 1*numPoints] = rlam3m5 + cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;

            index_tmp = 1*sz2 + i;
            for (k=5; k<nch; k++)
                for (j=0; j<5; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=5; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1 = -gam1 * (s1 - rlam3) * vv / c2 + (s2 * ny / c);
        cc2 = -gam1 * s2 * vv / c + (s1 - rlam3) * ny;

        index_tmp = 2*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = cc1 * uv + cc2 * nx;
        tau[index_tmp + 2*numPoints] = rlam3 + cc1 * vv + cc2 * ny;
        tau[index_tmp + 3*numPoints] = cc1 * wv + cc2 * nz;
        tau[index_tmp + 4*numPoints] = cc1 * h + cc2 * un;
        for (j=5; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = -gam1 * ((s1m1 - rlam3m1) * vv / c2 + (s1 - rlam3) * vvm1 / c2 - (s1 - rlam3) * vv * c2m1 / (c2 * c2)) +
                    s2m1 * ny / c - s2 * ny * cm1 / (c * c);
            cc1m2 = -gam1 * ((s1m2 - rlam3m2) * vv / c2 + (s1 - rlam3) * vvm2 / c2 - (s1 - rlam3) * vv * c2m2 / (c2 * c2)) +
                    s2m2 * ny / c - s2 * ny * cm2 / (c * c);
            cc1m3 = -gam1 * ((s1m3 - rlam3m3) * vv / c2 + (s1 - rlam3) * vvm3 / c2 - (s1 - rlam3) * vv * c2m3 / (c2 * c2)) +
                    s2m3 * ny / c - s2 * ny * cm3 / (c * c);
            cc1m4 = -gam1 * ((s1m4 - rlam3m4) * vv / c2 + (s1 - rlam3) * vvm4 / c2 - (s1 - rlam3) * vv * c2m4 / (c2 * c2)) +
                    s2m4 * ny / c - s2 * ny * cm4 / (c * c);
            cc1m5 = -gam1 * ((s1m5 - rlam3m5) * vv / c2 + (s1 - rlam3) * vvm5 / c2 - (s1 - rlam3) * vv * c2m5 / (c2 * c2)) +
                    s2m5 * ny / c - s2 * ny * cm5 / (c * c);

            cc2m1 = -gam1 * (s2m1 * vv / c + s2 * vvm1 / c - s2 * vv * cm1 / (c * c)) + (s1m1 - rlam3m1) * ny;
            cc2m2 = -gam1 * (s2m2 * vv / c + s2 * vvm2 / c - s2 * vv * cm2 / (c * c)) + (s1m2 - rlam3m2) * ny;
            cc2m3 = -gam1 * (s2m3 * vv / c + s2 * vvm3 / c - s2 * vv * cm3 / (c * c)) + (s1m3 - rlam3m3) * ny;
            cc2m4 = -gam1 * (s2m4 * vv / c + s2 * vvm4 / c - s2 * vv * cm4 / (c * c)) + (s1m4 - rlam3m4) * ny;
            cc2m5 = -gam1 * (s2m5 * vv / c + s2 * vvm5 / c - s2 * vv * cm5 / (c * c)) + (s1m5 - rlam3m5) * ny;

            index_tmp = (0*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m1 + cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m2 + cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m3 + cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m4 + cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            index_tmp += sz3; // (4*nch2 + 2*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m5;
            tau_uh[index_tmp + 1*numPoints] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            tau_uh[index_tmp + 2*numPoints] = rlam3m5 + cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;

            index_tmp = 2*sz2 + i;
            for (k=5; k<nch; k++)
                for (j=0; j<5; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=5; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1 = -gam1 * (s1 - rlam3) * wv / c2 + (s2 * nz / c);
        cc2 = -gam1 * s2 * wv / c + (s1 - rlam3) * nz;

        index_tmp = 3*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = cc1 * uv + cc2 * nx;
        tau[index_tmp + 2*numPoints] = cc1 * vv + cc2 * ny;
        tau[index_tmp + 3*numPoints] = rlam3 + cc1 * wv + cc2 * nz;
        tau[index_tmp + 4*numPoints] = cc1 * h + cc2 * un;
        for (j=5; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = -gam1 * ((s1m1 - rlam3m1) * wv / c2 + (s1 - rlam3) * wvm1 / c2 - (s1 - rlam3) * wv * c2m1 / (c2 * c2)) +
                    s2m1 * nz / c - s2 * nz * cm1 / (c * c);
            cc1m2 = -gam1 * ((s1m2 - rlam3m2) * wv / c2 + (s1 - rlam3) * wvm2 / c2 - (s1 - rlam3) * wv * c2m2 / (c2 * c2)) +
                    s2m2 * nz / c - s2 * nz * cm2 / (c * c);
            cc1m3 = -gam1 * ((s1m3 - rlam3m3) * wv / c2 + (s1 - rlam3) * wvm3 / c2 - (s1 - rlam3) * wv * c2m3 / (c2 * c2)) +
                    s2m3 * nz / c - s2 * nz * cm3 / (c * c);
            cc1m4 = -gam1 * ((s1m4 - rlam3m4) * wv / c2 + (s1 - rlam3) * wvm4 / c2 - (s1 - rlam3) * wv * c2m4 / (c2 * c2)) +
                    s2m4 * nz / c - s2 * nz * cm4 / (c * c);
            cc1m5 = -gam1 * ((s1m5 - rlam3m5) * wv / c2 + (s1 - rlam3) * wvm5 / c2 - (s1 - rlam3) * wv * c2m5 / (c2 * c2)) +
                    s2m5 * nz / c - s2 * nz * cm5 / (c * c);

            cc2m1 = -gam1 * (s2m1 * wv / c + s2 * wvm1 / c - s2 * wv * cm1 / (c * c)) + (s1m1 - rlam3m1) * nz;
            cc2m2 = -gam1 * (s2m2 * wv / c + s2 * wvm2 / c - s2 * wv * cm2 / (c * c)) + (s1m2 - rlam3m2) * nz;
            cc2m3 = -gam1 * (s2m3 * wv / c + s2 * wvm3 / c - s2 * wv * cm3 / (c * c)) + (s1m3 - rlam3m3) * nz;
            cc2m4 = -gam1 * (s2m4 * wv / c + s2 * wvm4 / c - s2 * wv * cm4 / (c * c)) + (s1m4 - rlam3m4) * nz;
            cc2m5 = -gam1 * (s2m5 * wv / c + s2 * wvm5 / c - s2 * wv * cm5 / (c * c)) + (s1m5 - rlam3m5) * nz;

            index_tmp = (0*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m1 + cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m2 + cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m3 + cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m4 + cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            index_tmp += sz3; // (4*nch2 + 3*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m5;
            tau_uh[index_tmp + 1*numPoints] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            tau_uh[index_tmp + 3*numPoints] = rlam3m5 + cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            tau_uh[index_tmp + 4*numPoints] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;

            index_tmp = 3*sz2 + i;
            for (k=5; k<nch; k++)
                for (j=0; j<5; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=5; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        cc1 = gam1 * (s1 - rlam3) / c2;
        cc2 = gam1 * s2 / c;

        index_tmp = 4*nch*numPoints+i;
        tau[index_tmp + 0*numPoints] = cc1;
        tau[index_tmp + 1*numPoints] = cc1 * uv + cc2 * nx;
        tau[index_tmp + 2*numPoints] = cc1 * vv + cc2 * ny;
        tau[index_tmp + 3*numPoints] = cc1 * wv + cc2 * nz;
        tau[index_tmp + 4*numPoints] = rlam3 + cc1 * h + cc2 * un;
        for (j=5; j<nch; j++)
            tau[index_tmp + j*numPoints] = 0.0;

        if (computeJacobian == 1) {
            cc1m1 = gam1 * ((s1m1 - rlam3m1) / c2 - (s1 - rlam3) * c2m1 / (c2 * c2));
            cc1m2 = gam1 * ((s1m2 - rlam3m2) / c2 - (s1 - rlam3) * c2m2 / (c2 * c2));
            cc1m3 = gam1 * ((s1m3 - rlam3m3) / c2 - (s1 - rlam3) * c2m3 / (c2 * c2));
            cc1m4 = gam1 * ((s1m4 - rlam3m4) / c2 - (s1 - rlam3) * c2m4 / (c2 * c2));
            cc1m5 = gam1 * ((s1m5 - rlam3m5) / c2 - (s1 - rlam3) * c2m5 / (c2 * c2));

            cc2m1 = gam1 * (s2m1 / c - s2 * cm1 / (c * c));
            cc2m2 = gam1 * (s2m2 / c - s2 * cm2 / (c * c));
            cc2m3 = gam1 * (s2m3 / c - s2 * cm3 / (c * c));
            cc2m4 = gam1 * (s2m4 / c - s2 * cm4 / (c * c));
            cc2m5 = gam1 * (s2m5 / c - s2 * cm5 / (c * c));

            index_tmp = (0*nch2 + 4*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m1;
            tau_uh[index_tmp + 1*numPoints] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            tau_uh[index_tmp + 4*numPoints] = rlam3m1 + cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            index_tmp += sz3; // (1*nch2 + 4*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m2;
            tau_uh[index_tmp + 1*numPoints] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            tau_uh[index_tmp + 4*numPoints] = rlam3m2 + cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            index_tmp += sz3; // (2*nch2 + 4*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m3;
            tau_uh[index_tmp + 1*numPoints] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            tau_uh[index_tmp + 4*numPoints] = rlam3m3 + cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            index_tmp += sz3; // (3*nch2 + 4*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m4;
            tau_uh[index_tmp + 1*numPoints] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            tau_uh[index_tmp + 4*numPoints] = rlam3m4 + cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            index_tmp += sz3; // (4*nch2 + 4*nch) * numPoints + i;
            tau_uh[index_tmp + 0*numPoints] = cc1m5;
            tau_uh[index_tmp + 1*numPoints] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            tau_uh[index_tmp + 2*numPoints] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            tau_uh[index_tmp + 3*numPoints] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            tau_uh[index_tmp + 4*numPoints] = rlam3m5 + cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;

            index_tmp = 4*sz2 + i;
            for (k=5; k<nch; k++)
                for (j=0; j<5; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;

            for (k=0; k<nch; k++)
                for (j=5; j<nch; j++)
                    tau_uh[index_tmp + k*sz3 + j*numPoints] = 0.0;
        }

        for (k=5; k<nch; k++)
            for (j=0; j<nch; j++)
                tau[k*sz2 + j*numPoints + i] = (k == j) ? stabTurbEq : 0.0;

        if (computeJacobian == 1) {
            for (m = 0; m < nch; m++)
                for (k = 5; k < nch; k++)
                    for (j = 0; j < nch; j++)
                        tau_uh[m*sz3 + k*sz2 + j*numPoints + i] = (k == j) ? DstabTurbEq[m] : 0.0;
        }
    }
    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch
    
    // Make zero the drivatives:
    for (i=0; i<numPoints*nch*nch*nch; i++)
        tau_uh[i] = 0.0;

    delete[] DstabTurbEq;
}


void getUpwindedStabilizationTensor(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    if (nd == 2)
        getUpwindedStabilizationTensor2d(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (nd == 3)
        getUpwindedStabilizationTensor3d(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else
        error("Number of dimensions not implemented.\n");
    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch
}


void getPropose3StabilizationTensor2d(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    // This functions works for compressible Euler, NS and RANS (any turbulence model).
    // If tauValue <= 0, the maximum eigenvalue is used to stabilize the turbulence model equation(s). Otherwise, tauValue is used.

//    getan2d(tau, tau_uh, uh, nl, param, 1, numPoints);

    double gam, gam1, epslm, tauValue, signun, ic, Dc2Reg_Dc2;
    double nx, ny;
    double r, ru, rv, rE;
    double r1, r1m1, r1m2, r1m3, r1m4;
    double uv, uvm1, uvm2, uvm3, uvm4;
    double vv, vvm1, vvm2, vvm3, vvm4;
    double Vgx, Vgy, Vgn;
    double E, Em1, Em2, Em3, Em4;
    double af, afm1, afm2, afm3, afm4;
    double p, pm1, pm2, pm3, pm4;
    double h, hm1, hm2, hm3, hm4;
    double c2, c2m1, c2m2, c2m3, c2m4;
    double c, cm1, cm2, cm3, cm4;
    double un, unm1, unm2, unm3, unm4;
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4;
    double rlam1, rlam1m1, rlam1m2, rlam1m3, rlam1m4;
    double rlam2, rlam2m1, rlam2m2, rlam2m3, rlam2m4;
    double lamMax, lamMaxm1, lamMaxm2, lamMaxm3, lamMaxm4;

    double eps = 1.0e-8;        // eps > 0 ensures that tau is non-zero on at least one side of all faces even with finite-precision arithmetic
    double pi = 3.141592653589793;
    double b = 100.0;
    
    gam   = param[0];
    gam1 = gam - 1.0;
    epslm = param[1];
    tauValue = param[5];

    double stabTurbEq;
    double *DstabTurbEq  = new double [nch];

    int nch2 = nch * nch, sz2 = nch * numPoints, sz3 = sz2 * nch;
    int i, j, k, m, index_tmp;

    Int ALEflag = app.ALEflag;

    for (i=0; i<numPoints; i++) {
        r   = uh[0*numPoints+i]; // 1.0000;
        ru  = uh[1*numPoints+i]; // 0.9994;
        rv  = uh[2*numPoints+i]; // 0.0349;
        rE  = uh[3*numPoints+i]; // 45.1429;

        nx  = nl[0*numPoints+i];
        ny  = nl[1*numPoints+i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
        }

        r1   = 1.0/r;
        uv   = ru*r1;
        vv   = rv*r1;
        E    = rE*r1;
        af   = 0.5*(uv*uv+vv*vv);
        p    = gam1*(rE -r*af);
        h    = E   + p*r1;
        c2   = gam* p*r1;
        Dc2Reg_Dc2 = atan(b*c2)/pi + 1.0/2.0 + (b*c2) / ((1+b*c2*b*c2) * pi);
        c2 = max ( c2*(atan(b*c2)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 , 1.0e-6);
        c    = sqrt(c2);
        un   = uv*nx   + vv*ny;

        r1m1 = -1.0/(r*r);
        r1m2 = 0.0;
        r1m3 = 0.0;
        r1m4 = 0.0;

        if (computeJacobian == 1) {
            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;
            uvm3 = ru * r1m3;
            uvm4 = ru * r1m4;

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;
            vvm3 = r1 + rv * r1m3;
            vvm4 = rv * r1m4;

            Em1 = rE * r1m1;
            Em2 = rE * r1m2;
            Em3 = rE * r1m3;
            Em4 = r1 + rE * r1m4;

            afm1 = uv * uvm1 + vv * vvm1;
            afm2 = uv * uvm2 + vv * vvm2;
            afm3 = uv * uvm3 + vv * vvm3;
            afm4 = uv * uvm4 + vv * vvm4;

            pm1 = gam1 * (-af - r * afm1);
            pm2 = gam1 * (-r * afm2);
            pm3 = gam1 * (-r * afm3);
            pm4 = gam1 * (1.0 - r * afm4);

            hm1 = Em1 + pm1 * r1 + p * r1m1;
            hm2 = Em2 + pm2 * r1 + p * r1m2;
            hm3 = Em3 + pm3 * r1 + p * r1m3;
            hm4 = Em4 + pm4 * r1 + p * r1m4;

            c2m1 = gam * (pm1 * r1 + p * r1m1);
            c2m2 = gam * (pm2 * r1 + p * r1m2);
            c2m3 = gam * (pm3 * r1 + p * r1m3);
            c2m4 = gam * (pm4 * r1 + p * r1m4);
            
            c2m1 *= Dc2Reg_Dc2;
            c2m2 *= Dc2Reg_Dc2;
            c2m3 *= Dc2Reg_Dc2;
            c2m4 *= Dc2Reg_Dc2;

            cm1 = 0.5 * c2m1 / c;
            cm2 = 0.5 * c2m2 / c;
            cm3 = 0.5 * c2m3 / c;
            cm4 = 0.5 * c2m4 / c;

            unm1 = uvm1 * nx + vvm1 * ny;
            unm2 = uvm2 * nx + vvm2 * ny;
            unm3 = uvm3 * nx + vvm3 * ny;
            unm4 = uvm4 * nx + vvm4 * ny;
        }
        
        Vgn = Vgx * nx + Vgy * ny;

        rlam1 = (un - Vgn) + fabs(un - Vgn);
        signun = (un - Vgn >= 0) ? 1.0 : -1.0;
        if (computeJacobian == 1) {
            rlam1m1 = (1 + signun) * unm1;
            rlam1m2 = (1 + signun) * unm2;
            rlam1m3 = (1 + signun) * unm3;
            rlam1m4 = (1 + signun) * unm4;
        }

        rlam2 = fabs(un - Vgn) + c;
        signun = (un - Vgn >= 0) ? 1.0 : -1.0;
        if (computeJacobian == 1) {
            rlam2m1 = signun * unm1 + cm1;
            rlam2m2 = signun * unm2 + cm2;
            rlam2m3 = signun * unm3 + cm3;
            rlam2m4 = signun * unm4 + cm4;
        }

        if (un - Vgn >= 0) {
            lamMax = (un-Vgn) + c;
            if (computeJacobian == 1) {
                lamMaxm1 = unm1 + cm1;
                lamMaxm2 = unm2 + cm2;
                lamMaxm3 = unm3 + cm3;
                lamMaxm4 = unm4 + cm4;
            }
        }
        else {
            lamMax = - (un-Vgn) + c;
            if (computeJacobian == 1) {
                lamMaxm1 = -unm1 + cm1;
                lamMaxm2 = -unm2 + cm2;
                lamMaxm3 = -unm3 + cm3;
                lamMaxm4 = -unm4 + cm4;
            }
        }

        if (tauValue > 0) {
            stabTurbEq = tauValue;
            if (computeJacobian == 1) {
                for (j = 0; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }
        else {
            stabTurbEq = lamMax;
            if (computeJacobian == 1) {
                DstabTurbEq[0] = lamMaxm1;
                DstabTurbEq[1] = lamMaxm2;
                DstabTurbEq[2] = lamMaxm3;
                DstabTurbEq[3] = lamMaxm4;
                for (j = 4; j < nch; j++)
                    DstabTurbEq[j] = 0.0;
            }
        }

        for (j=0; j<nch*nch; j++)
            tau[i + j*numPoints] = 0.0;
        tau[i + 0*numPoints]  = rlam1;
        tau[i + 5*numPoints]  = rlam1;
        tau[i + 10*numPoints]  = rlam1;
        tau[i + 15*numPoints]  = rlam2;

        if (computeJacobian == 1) {
            for (j=0; j<nch*nch*nch; j++)
                tau_uh[i + j*numPoints] = 0.0;

            tau_uh[i + 0*numPoints + 0*sz3] = rlam1m1;
            tau_uh[i + 0*numPoints + 1*sz3] = rlam1m2;
            tau_uh[i + 0*numPoints + 2*sz3] = rlam1m3;
            tau_uh[i + 0*numPoints + 3*sz3] = rlam1m4;

            tau_uh[i + 5*numPoints + 0*sz3] = rlam1m1;
            tau_uh[i + 5*numPoints + 1*sz3] = rlam1m2;
            tau_uh[i + 5*numPoints + 2*sz3] = rlam1m3;
            tau_uh[i + 5*numPoints + 3*sz3] = rlam1m4;

            tau_uh[i + 10*numPoints + 0*sz3] = rlam1m1;
            tau_uh[i + 10*numPoints + 1*sz3] = rlam1m2;
            tau_uh[i + 10*numPoints + 2*sz3] = rlam1m3;
            tau_uh[i + 10*numPoints + 3*sz3] = rlam1m4;

            tau_uh[i + 15*numPoints + 0*sz3] = rlam2m1;
            tau_uh[i + 15*numPoints + 1*sz3] = rlam2m2;
            tau_uh[i + 15*numPoints + 2*sz3] = rlam2m3;
            tau_uh[i + 15*numPoints + 3*sz3] = rlam2m4;
        }

        for (j=4; j<nch; j++)
            tau[j*sz2 + j*numPoints + i] = stabTurbEq;

        if (computeJacobian == 1) {
            for (m = 0; m < nch; m++)
                for (j = 4; j < nch; j++)
                    tau_uh[m*sz3 + j*sz2 + j*numPoints + i] = DstabTurbEq[m];
        }
    }
    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch

   for (i=0; i<numPoints*nch*nch*nch; i++)
       tau_uh[i] = 0.0;

    delete[] DstabTurbEq;
}

void getPropose3StabilizationTensor(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    if (nd == 2)
        getPropose3StabilizationTensor2d(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (nd == 3) {
        error("No implemented.\n");
//        getPropose3StabilizationTensor3d(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    }
    else
        error("Number of dimensions not implemented.\n");
    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch
}

void getConvectiveStabilizationTensor(double* tau, double* tau_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, Int convStabMethod, int numPoints, int nch, int nd, int computeJacobian)
{
    if (convStabMethod == 0)
        getConstantStabilizationTensor(tau, tau_uh, param, numPoints, nch, computeJacobian);
    else if (convStabMethod == 1)
        getLaxFriedrichStabilizationTensor(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (convStabMethod == 2)
        getRoeStabilizationTensor(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (convStabMethod == 3)
        getProposedStabilizationTensor(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (convStabMethod == 4)
        getUpwindedStabilizationTensor(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (convStabMethod == 5)
        getPropose3StabilizationTensor(tau, tau_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else
        error("Stabilization method for convective operator not implemented.\n");
}

void getConstantDiffStabilizationTensor(double* tauDiff, double* tauDiff_uh, double* uh, double* pg, double* nl, 
        appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    if (nch != 2+nd)
        error("getTensorialDiffStabilizationTensor not implemented for RANS equations.\n");
    
    int i, j, k;
    int sz2 = nch * numPoints, sz3 = sz2 * nch, sz4 = sz3 * nch;
    
    double Re = param[2];
    double kappa = 1.0 / Re;             // kappa = 1 / Re
    
    for (k=0; k<nch; k++)
        for (j=0; j<nch; j++) {
            for (i = 0; i < numPoints; i++) {
                double l = 1.0 / sqrt(Re); //pg[?? * numPoints + i];
                tauDiff[i + j * numPoints + k * sz2] = (j==k) ? (kappa / l) : 0.0;
            }
        }
    // tauDiff: numPoints / nch / nch
    
    if (computeJacobian == 1) {
        for (i=0; i<sz4; i++)
            tauDiff_uh[i] = 0.0;
    }
    // tauDiff_uh: numPoints / nch / nch / nch
}

void getTensorialDiffStabilizationTensor2d(double* tauDiff, double* tauDiff_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
//     error("getTensorialDiffStabilizationTensor2d not validated yet.\n");
    
    if (nch != 2+nd)
        error("getTensorialDiffStabilizationTensor2d not implemented for RANS equations.\n");
    
    int i, sz, nch2 = nch * nch, sz2 = nch * numPoints, sz3 = sz2 * nch, sz4 = sz3 * nch;
    
    double nx, ny;
    double r, ru, rv, rE;
    double Vgx, Vgy, Vgn;
    double r1, r1m1, r1m2, r1m3, r1m4;
    double uv, uvm1, uvm2, uvm3, uvm4;
    double vv, vvm1, vvm2, vvm3, vvm4;
    double E, Em1, Em2, Em3, Em4;
    double af, afm1, afm2, afm3, afm4;
    double un, unm1, unm2, unm3, unm4;
    double uvuv, uvuvm1, uvuvm2, uvuvm3, uvuvm4;
    double uvvv, uvvvm1, uvvvm2, uvvvm3, uvvvm4;
    double vvuv, vvuvm1, vvuvm2, vvuvm3, vvuvm4;
    double vvvv, vvvvm1, vvvvm2, vvvvm3, vvvvm4;
    double uvuvr1m1, uvuvr1m2, uvuvr1m3, uvuvr1m4;
    double uvvvr1m1, uvvvr1m2, uvvvr1m3, uvvvr1m4;
    double vvuvr1m1, vvuvr1m2, vvuvr1m3, vvuvr1m4;
    double vvvvr1m1, vvvvr1m2, vvvvr1m3, vvvvr1m4;
    double uvr1m1, uvr1m2, uvr1m3, uvr1m4;
    double vvr1m1, vvr1m2, vvr1m3, vvr1m4;
    double afr1m1, afr1m2, afr1m3, afr1m4;
    double unr1m1, unr1m2, unr1m3, unr1m4;
    double Er1m1, Er1m2, Er1m3, Er1m4;
    
    double gam = param[0];
    double gam1 = gam - 1.0;
    double Re = param[2];
    double Pr = param[3];
    double pi = 3.141592653589793;
    double b = 100.0;
    
    double l_ref;

    Int ALEflag = app.ALEflag;

    for (i = 0; i < sz3; i++)
        tauDiff[i] = 0.0;
    for (i = 0; i < sz4; i++)
        tauDiff_uh[i] = 0.0;
    
    for (i = 0; i < numPoints; i++) {
        nx = nl[0 * numPoints + i];
        ny = nl[1 * numPoints + i];

        r = uh[0 * numPoints + i];
        ru = uh[1 * numPoints + i];
        rv = uh[2 * numPoints + i];
        rE = uh[3 * numPoints + i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
        }
        Vgn = Vgx * nx + Vgy * ny;
        
        r1 = 1.0 / r;
        uv = ru * r1;
        vv = rv * r1;
        E = rE * r1;
        af = 0.5 * (uv * uv + vv * vv);
        un = uv * nx + vv * ny;
        
        l_ref = 1.0 / sqrt(Re); //pg[?? * numPoints + i];
        
        tauDiff[i + 1 * numPoints + 0 * sz2] = - (1.0 / Re) * ( uv * r1 + nx * un / (3*r) ) / l_ref;
        tauDiff[i + 2 * numPoints + 0 * sz2] = - (1.0 / Re) * ( vv * r1 + ny * un / (3*r) ) / l_ref;

        tauDiff[i + 1 * numPoints + 1 * sz2] = (1.0 / Re) * ( r1 + nx * nx / (3*r) ) / l_ref;
        tauDiff[i + 2 * numPoints + 1 * sz2] = (1.0 / Re) * (      nx * ny / (3*r) ) / l_ref;

        tauDiff[i + 1 * numPoints + 2 * sz2] = (1.0 / Re) * (      ny * nx / (3*r) ) / l_ref;
        tauDiff[i + 2 * numPoints + 2 * sz2] = (1.0 / Re) * ( r1 + ny * ny / (3*r) ) / l_ref;

        tauDiff[i + (1+nd) * numPoints + 0 * sz2] = - (1.0 / Re) * ( 2.0*af * r1 + (nx*uv*uv*nx + nx*uv*vv*ny + ny*vv*uv*nx + ny*vv*vv*ny) / (3*r) ) / l_ref
                                                    + ( - (gam/(Re*Pr)) * E * r1 + (gam/(Re*Pr)) * 2.0*af * r1 ) / l_ref;
        
        tauDiff[i + (1+nd) * numPoints + 1 * sz2] = (1.0 / Re) * ( uv * r1 + nx*un / (3*r) ) / l_ref - ((gam/(Re*Pr)) * uv * r1) / l_ref;
        tauDiff[i + (1+nd) * numPoints + 2 * sz2] = (1.0 / Re) * ( vv * r1 + ny*un / (3*r) ) / l_ref - ((gam/(Re*Pr)) * vv * r1) / l_ref;

        tauDiff[i + (1+nd) * numPoints + (1+nd) * sz2] = ((gam/(Re*Pr)) * r1) / l_ref;

        if (computeJacobian == 1) {
            r1m1 = -1.0 / (r * r);
            r1m2 = 0.0;
            r1m3 = 0.0;
            r1m4 = 0.0;

            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;//r1
            uvm3 = ru * r1m3;//0.0
            uvm4 = ru * r1m4;//0.0

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;//0.0
            vvm3 = r1 + rv * r1m3;//r1
            vvm4 = rv * r1m4;//0.0
            
            uvuvm1 = 2.0*uv*uvm1;
            uvuvm2 = 2.0*uv*uvm2;
            uvuvm3 = 2.0*uv*uvm3;//0.0
            uvuvm4 = 2.0*uv*uvm4;//0.0
            
            uvvvm1 = uv*vvm1 + uvm1*vv;
            uvvvm2 = uv*vvm2 + uvm2*vv;
            uvvvm3 = uv*vvm3 + uvm3*vv;
            uvvvm4 = uv*vvm4 + uvm4*vv;//0.0
            
            vvuvm1 = vv*uvm1 + vvm1*uv;
            vvuvm2 = vv*uvm2 + vvm2*uv;
            vvuvm3 = vv*uvm3 + vvm3*uv;
            vvuvm4 = vv*uvm4 + vvm4*uv;//0.0
            
            vvvvm1 = 2.0*vv*vvm1;
            vvvvm2 = 2.0*vv*vvm2;//0.0
            vvvvm3 = 2.0*vv*vvm3;
            vvvvm4 = 2.0*vv*vvm4;//0.0
            
            Em1 = rE * r1m1;
            Em2 = rE * r1m2;//0.0
            Em3 = rE * r1m3;//0.0
            Em4 = rE * r1m4;

            afm1 = uv * uvm1 + vv * vvm1;
            afm2 = uv * uvm2 + vv * vvm2;
            afm3 = uv * uvm3 + vv * vvm3;
            afm4 = uv * uvm4 + vv * vvm4;//0.0

            unm1 = uvm1 * nx + vvm1 * ny;
            unm2 = uvm2 * nx + vvm2 * ny;
            unm3 = uvm3 * nx + vvm3 * ny;
            unm4 = uvm4 * nx + vvm4 * ny;//0.0
            
            Er1m1 = Em1*r1+E*r1m1;
            Er1m2 = Em2*r1+E*r1m2;
            Er1m3 = Em3*r1+E*r1m3;
            Er1m4 = Em4*r1+E*r1m4;
            
            uvr1m1 = uvm1*r1+uv*r1m1;
            uvr1m2 = uvm2*r1+uv*r1m2;
            uvr1m3 = uvm3*r1+uv*r1m3;//0.0
            uvr1m4 = uvm4*r1+uv*r1m4;//0.0
            
            vvr1m1 = vvm1*r1+vv*r1m1;
            vvr1m2 = vvm2*r1+vv*r1m2;//0.0
            vvr1m3 = vvm3*r1+vv*r1m3;
            vvr1m4 = vvm4*r1+vv*r1m4;//0.0
            
            uvuvr1m1 = uvuvm1*r1 + uvuv*r1m1;
            uvuvr1m2 = uvuvm2*r1 + uvuv*r1m2;
            uvuvr1m3 = uvuvm3*r1 + uvuv*r1m3;
            uvuvr1m4 = uvuvm4*r1 + uvuv*r1m4;
            
            uvvvr1m1 = uvvvm1*r1 + uvvv*r1m1;
            uvvvr1m2 = uvvvm2*r1 + uvvv*r1m2;
            uvvvr1m3 = uvvvm3*r1 + uvvv*r1m3;
            uvvvr1m4 = uvvvm4*r1 + uvvv*r1m4;
            
            vvuvr1m1 = vvuvm1*r1 + vvuv*r1m1;
            vvuvr1m2 = vvuvm2*r1 + vvuv*r1m2;
            vvuvr1m3 = vvuvm3*r1 + vvuv*r1m3;
            vvuvr1m4 = vvuvm4*r1 + vvuv*r1m4;
            
            vvvvr1m1 = vvvvm1*r1 + vvvv*r1m1;
            vvvvr1m2 = vvvvm2*r1 + vvvv*r1m2;
            vvvvr1m3 = vvvvm3*r1 + vvvv*r1m3;
            vvvvr1m4 = vvvvm4*r1 + vvvv*r1m4;
            
            unr1m1 = unm1*r1+un*r1m1;
            unr1m2 = unm2*r1+un*r1m2;
            unr1m3 = unm3*r1+un*r1m3;
            unr1m4 = unm4*r1+un*r1m4;//0.0
            
            afr1m1 = afm1*r1+af*r1m1;
            afr1m2 = afm2*r1+af*r1m2;
            afr1m3 = afm3*r1+af*r1m3;
            afr1m4 = afm4*r1+af*r1m4;//0.0
            
            tauDiff_uh[i + 1 * numPoints + 0 * sz2 + 0 * sz3] = - (1.0 / Re) * ( uvr1m1 + nx * unr1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 0 * sz2 + 1 * sz3] = - (1.0 / Re) * ( uvr1m2 + nx * unr1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 0 * sz2 + 2 * sz3] = - (1.0 / Re) * ( uvr1m3 + nx * unr1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 0 * sz2 + 3 * sz3] = - (1.0 / Re) * ( uvr1m4 + nx * unr1m4 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 2 * numPoints + 0 * sz2 + 0 * sz3] = - (1.0 / Re) * ( vvr1m1 + ny * unr1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 0 * sz2 + 1 * sz3] = - (1.0 / Re) * ( vvr1m2 + ny * unr1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 0 * sz2 + 2 * sz3] = - (1.0 / Re) * ( vvr1m3 + ny * unr1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 0 * sz2 + 3 * sz3] = - (1.0 / Re) * ( vvr1m4 + ny * unr1m4 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 1 * numPoints + 1 * sz2 + 0 * sz3] = (1.0 / Re) * ( r1m1 + nx * nx * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 1 * sz2 + 1 * sz3] = (1.0 / Re) * ( r1m2 + nx * nx * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 1 * sz2 + 2 * sz3] = (1.0 / Re) * ( r1m3 + nx * nx * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 1 * sz2 + 3 * sz3] = (1.0 / Re) * ( r1m4 + nx * nx * r1m4 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 2 * numPoints + 1 * sz2 + 0 * sz3] = (1.0 / Re) * (           nx * ny * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 1 * sz2 + 1 * sz3] = (1.0 / Re) * (           nx * ny * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 1 * sz2 + 2 * sz3] = (1.0 / Re) * (           nx * ny * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 1 * sz2 + 3 * sz3] = (1.0 / Re) * (           nx * ny * r1m4 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 1 * numPoints + 2 * sz2 + 0 * sz3] = (1.0 / Re) * (           ny * nx * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 2 * sz2 + 1 * sz3] = (1.0 / Re) * (           ny * nx * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 2 * sz2 + 2 * sz3] = (1.0 / Re) * (           ny * nx * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 2 * sz2 + 3 * sz3] = (1.0 / Re) * (           ny * nx * r1m4 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 2 * numPoints + 2 * sz2 + 0 * sz3] = (1.0 / Re) * ( r1m1 + ny * ny * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 2 * sz2 + 1 * sz3] = (1.0 / Re) * ( r1m2 + ny * ny * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 2 * sz2 + 2 * sz3] = (1.0 / Re) * ( r1m3 + ny * ny * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 2 * sz2 + 3 * sz3] = (1.0 / Re) * ( r1m4 + ny * ny * r1m4 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + (1+nd) * numPoints + 0 * sz2 + 0 * sz3] = - (1.0 / Re) * ( 2.0*afr1m1 + (nx*uvuvr1m1*nx + nx*uvvvr1m1*ny + ny*vvuvr1m1*nx + ny*vvvvr1m1*ny) / 3.0 ) / l_ref
                                                                    + ( - (gam/(Re*Pr)) * Er1m1 + (gam/(Re*Pr)) * 2.0*afr1m1 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 0 * sz2 + 1 * sz3] = - (1.0 / Re) * ( 2.0*afr1m2 + (nx*uvuvr1m2*nx + nx*uvvvr1m2*ny + ny*vvuvr1m2*nx + ny*vvvvr1m2*ny) / 3.0 ) / l_ref
                                                                    + ( - (gam/(Re*Pr)) * Er1m2 + (gam/(Re*Pr)) * 2.0*afr1m2 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 0 * sz2 + 2 * sz3] = - (1.0 / Re) * ( 2.0*afr1m3 + (nx*uvuvr1m3*nx + nx*uvvvr1m3*ny + ny*vvuvr1m3*nx + ny*vvvvr1m3*ny) / 3.0 ) / l_ref
                                                                    + ( - (gam/(Re*Pr)) * Er1m3 + (gam/(Re*Pr)) * 2.0*afr1m3 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 0 * sz2 + 3 * sz3] = - (1.0 / Re) * ( 2.0*afr1m4 + (nx*uvuvr1m4*nx + nx*uvvvr1m4*ny + ny*vvuvr1m4*nx + ny*vvvvr1m4*ny) / 3.0 ) / l_ref
                                                                    + ( - (gam/(Re*Pr)) * Er1m4 + (gam/(Re*Pr)) * 2.0*afr1m4 ) / l_ref;
            
            tauDiff_uh[i + (1+nd) * numPoints + 1 * sz2 + 0 * sz3] = (1.0 / Re) * ( uvr1m1 + nx*unr1m1 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * uvr1m1 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 1 * sz2 + 1 * sz3] = (1.0 / Re) * ( uvr1m2 + nx*unr1m2 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * uvr1m2 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 1 * sz2 + 2 * sz3] = (1.0 / Re) * ( uvr1m3 + nx*unr1m3 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * uvr1m3 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 1 * sz2 + 3 * sz3] = (1.0 / Re) * ( uvr1m4 + nx*unr1m4 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * uvr1m4 ) / l_ref;
            
            tauDiff_uh[i + (1+nd) * numPoints + 2 * sz2 + 0 * sz3] = (1.0 / Re) * ( vvr1m1 + ny*unr1m1 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * vvr1m1 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 2 * sz2 + 1 * sz3] = (1.0 / Re) * ( vvr1m2 + ny*unr1m2 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * vvr1m2 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 2 * sz2 + 2 * sz3] = (1.0 / Re) * ( vvr1m3 + ny*unr1m3 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * vvr1m3 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 2 * sz2 + 3 * sz3] = (1.0 / Re) * ( vvr1m4 + ny*unr1m4 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * vvr1m4 ) / l_ref;
            
            tauDiff_uh[i + (1+nd) * numPoints + (1+nd) * sz2 + 0 * sz3] = ((gam/(Re*Pr)) * r1m1 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + (1+nd) * sz2 + 1 * sz3] = ((gam/(Re*Pr)) * r1m2 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + (1+nd) * sz2 + 2 * sz3] = ((gam/(Re*Pr)) * r1m3 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + (1+nd) * sz2 + 3 * sz3] = ((gam/(Re*Pr)) * r1m4 ) / l_ref;
        }
    }
    // tauDiff: numPoints / nch / nch
    // tauDiff_uh: numPoints / nch / nch / nch
    
//    for (i=0; i<sz4; i++)
//        tauDiff_uh[i] = 0.0;
}

void getTensorialDiffStabilizationTensor3d(double* tauDiff, double* tauDiff_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    // Viscous diffusivity tensor:
    // K_{0,l}^{j,s} = 0, l = 0, ..., 1+nd
    // K_{1+i,0}^{j,s} = - (mu/Re) * (u_i d_{s,j} / r + u_j d_{s,i} / r - (2/3) u_s d_{i,j} / r ) , i=0, ..., nd-1
    // K_{1+i,1+l}^{j,s} = (mu/Re) * (d_{l,i} d_{s,j} / r + d_{l,j} d_{s,i} / r - (2/3) d_{i,j} d_{l,s} / r ) , i,l=0, ..., nd-1
    // K_{1+i,1+nd}^{j,s} = 0, i=0, ..., nd-1
    // K_{1+nd,0}^{j,s} = - (mu/Re) * ( u_s u_j / 3r + ||u||^2 d_{s,j} / r ) - (mu*gam/Re*Pr) * E d_{j,s} / r + (mu*gam/Re*Pr) ||v||^2 d_{j,s} / r
    // K_{1+nd,1+l}^{j,s} = (mu/Re) * ( u_s d_{l,j} / r + u_l d_{s,j} / r - (2/3) u_j d_{l,s} / r ) - (mu*gam/Re*Pr) * u_l d_{j,s} / r , l=0, ..., nd-1
    // K_{1+nd,1+nd}^{j,s} = (mu*gam/Re*Pr) * d_{j,s} / r
    // Also, for RANS-SA:
    // K_{2+nd,0}^{j,s} = - nu^2 (1+Psi) d_{j,s} / sigma
    // K_{2+nd,l}^{j,s} = 0, l=1,...,1+nd
    // K_{2+nd,2+nd}^{j,s} = nu (1+Psi) d_{j,s} / sigma
    
    // Stabilization tensor for diffusive operator:
    // l_ref tauDiff_{0,l} = 0, l = 0, ..., 1+nd
    // l_ref tauDiff_{1+i,0} = - (mu/Re) * (u_i / r + v_n n_i / 3r ) , i=0, ..., nd-1
    // l_ref tauDiff_{1+i,1+l} = (mu/Re) * (d_{i,l} / r + n_i n_l / 3r ) , i,l=0, ..., nd-1
    // l_ref tauDiff_{1+i,1+nd} = 0, i=0, ..., nd-1
    // l_ref tauDiff_{1+nd,0} = - (mu/Re) * ( n_s u_s u_j n_j / 3r + ||v||^2 / r ) - (mu*gam/Re*Pr) * E / r + (mu*gam/Re*Pr) ||v||^2 / r
    // l_ref tauDiff_{1+nd,1+l} = (mu/Re) * ( v_n n_l / 3r + u_l / r ) - (mu*gam/Re*Pr) * u_l / r , l=0, ..., nd-1
    // l_ref tauDiff_{1+nd,1+nd} = (mu*gam/Re*Pr) / r
    // Also, for RANS-SA:
    // l_ref tauDiff_{2+nd,0} = - nu^2 (1+Psi) / sigma
    // l_ref tauDiff_{2+nd,l} = 0, l=1,...,1+nd
    // l_ref tauDiff_{2+nd,2+nd} = nu (1+Psi) / sigma
    
    // Note: Here, mu denotes the nondimensional dynamic viscosity (i.e. mu / mu_ref = mu / mu_{inf})
    
    error("getTensorialDiffStabilizationTensor3d not validated yet.\n");
    
    if (nch != 2+nd)
        error("getTensorialDiffStabilizationTensor3d not implemented for RANS equations.\n");
    
    int i, sz, nch2 = nch * nch, sz2 = nch * numPoints, sz3 = sz2 * nch, sz4 = sz3 * nch;
    
    double nx, ny, nz;
    double r, ru, rv, rw, rE;
    double Vgx, Vgy, Vgz, Vgn;
    double r1, r1m1, r1m2, r1m3, r1m4, r1m5;
    double uv, uvm1, uvm2, uvm3, uvm4, uvm5;
    double vv, vvm1, vvm2, vvm3, vvm4, vvm5;
    double wv, wvm1, wvm2, wvm3, wvm4, wvm5;
    double E, Em1, Em2, Em3, Em4, Em5;
    double af, afm1, afm2, afm3, afm4, afm5;
    double un, unm1, unm2, unm3, unm4, unm5;
    double uvuv, uvuvm1, uvuvm2, uvuvm3, uvuvm4, uvuvm5;
    double uvvv, uvvvm1, uvvvm2, uvvvm3, uvvvm4, uvvvm5;
    double uvwv, uvwvm1, uvwvm2, uvwvm3, uvwvm4, uvwvm5;
    double vvuv, vvuvm1, vvuvm2, vvuvm3, vvuvm4, vvuvm5;
    double vvvv, vvvvm1, vvvvm2, vvvvm3, vvvvm4, vvvvm5;
    double vvwv, vvwvm1, vvwvm2, vvwvm3, vvwvm4, vvwvm5;
    double wvuv, wvuvm1, wvuvm2, wvuvm3, wvuvm4, wvuvm5;
    double wvvv, wvvvm1, wvvvm2, wvvvm3, wvvvm4, wvvvm5;
    double wvwv, wvwvm1, wvwvm2, wvwvm3, wvwvm4, wvwvm5;
    double uvuvr1m1, uvuvr1m2, uvuvr1m3, uvuvr1m4, uvuvr1m5;
    double uvvvr1m1, uvvvr1m2, uvvvr1m3, uvvvr1m4, uvvvr1m5;
    double uvwvr1m1, uvwvr1m2, uvwvr1m3, uvwvr1m4, uvwvr1m5;
    double vvuvr1m1, vvuvr1m2, vvuvr1m3, vvuvr1m4, vvuvr1m5;
    double vvvvr1m1, vvvvr1m2, vvvvr1m3, vvvvr1m4, vvvvr1m5;
    double vvwvr1m1, vvwvr1m2, vvwvr1m3, vvwvr1m4, vvwvr1m5;
    double wvuvr1m1, wvuvr1m2, wvuvr1m3, wvuvr1m4, wvuvr1m5;
    double wvvvr1m1, wvvvr1m2, wvvvr1m3, wvvvr1m4, wvvvr1m5;
    double wvwvr1m1, wvwvr1m2, wvwvr1m3, wvwvr1m4, wvwvr1m5;
    double uvr1m1, uvr1m2, uvr1m3, uvr1m4, uvr1m5;
    double vvr1m1, vvr1m2, vvr1m3, vvr1m4, vvr1m5;
    double wvr1m1, wvr1m2, wvr1m3, wvr1m4, wvr1m5;
    double afr1m1, afr1m2, afr1m3, afr1m4, afr1m5;
    double unr1m1, unr1m2, unr1m3, unr1m4, unr1m5;
    double Er1m1, Er1m2, Er1m3, Er1m4, Er1m5;
    
    double gam = param[0];
    double gam1 = gam - 1.0;
    double Re = param[2];
    double Pr = param[3];
    double pi = 3.141592653589793;
    double b = 100.0;
    
    double l_ref;

    Int ALEflag = app.ALEflag;

    for (i = 0; i < sz3; i++)
        tauDiff[i] = 0.0;
    for (i = 0; i < sz4; i++)
        tauDiff_uh[i] = 0.0;
    
    for (i = 0; i < numPoints; i++) {
        nx = nl[0 * numPoints + i];
        ny = nl[1 * numPoints + i];
        nz = nl[2 * numPoints + i];

        r = uh[0 * numPoints + i];
        ru = uh[1 * numPoints + i];
        rv = uh[2 * numPoints + i];
        rw = uh[3 * numPoints + i];
        rE = uh[4 * numPoints + i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
            Vgz = pg[(2 * nd + 2) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
            Vgz = 0.0;
        }
        Vgn = Vgx * nx + Vgy * ny + Vgz * nz;
        
        r1 = 1.0 / r;
        uv = ru * r1;
        vv = rv * r1;
        wv = rw * r1;
        E = rE * r1;
        af = 0.5 * (uv * uv + vv * vv + wv * wv);
        un = uv * nx + vv * ny + wv * nz;
        
        l_ref = 1.0 / sqrt(Re); //pg[?? * numPoints + i];
        
        tauDiff[i + 1 * numPoints + 0 * sz2] = - (1.0 / Re) * ( uv * r1 + nx * un / (3*r) ) / l_ref;
        tauDiff[i + 2 * numPoints + 0 * sz2] = - (1.0 / Re) * ( vv * r1 + ny * un / (3*r) ) / l_ref;
        tauDiff[i + 3 * numPoints + 0 * sz2] = - (1.0 / Re) * ( wv * r1 + nz * un / (3*r) ) / l_ref;
        
        tauDiff[i + 1 * numPoints + 1 * sz2] = (1.0 / Re) * ( r1 + nx * nx / (3*r) ) / l_ref;
        tauDiff[i + 2 * numPoints + 1 * sz2] = (1.0 / Re) * (      nx * ny / (3*r) ) / l_ref;
        tauDiff[i + 3 * numPoints + 1 * sz2] = (1.0 / Re) * (      nx * nz / (3*r) ) / l_ref;
        
        tauDiff[i + 1 * numPoints + 2 * sz2] = (1.0 / Re) * (      ny * nx / (3*r) ) / l_ref;
        tauDiff[i + 2 * numPoints + 2 * sz2] = (1.0 / Re) * ( r1 + ny * ny / (3*r) ) / l_ref;
        tauDiff[i + 3 * numPoints + 2 * sz2] = (1.0 / Re) * (      ny * nz / (3*r) ) / l_ref;
        
        tauDiff[i + 1 * numPoints + 3 * sz2] = (1.0 / Re) * (      nz * nx / (3*r) ) / l_ref;
        tauDiff[i + 2 * numPoints + 3 * sz2] = (1.0 / Re) * (      nz * ny / (3*r) ) / l_ref;
        tauDiff[i + 3 * numPoints + 3 * sz2] = (1.0 / Re) * ( r1 + nz * nz / (3*r) ) / l_ref;
        
        tauDiff[i + (1+nd) * numPoints + 0 * sz2] = - (1.0 / Re) * ( 2.0*af * r1 + (nx*uv*uv*nx + nx*uv*vv*ny + nx*uv*wv*nz + ny*vv*uv*nx + ny*vv*vv*ny + ny*vv*wv*nz + nz*wv*uv*nx + nz*wv*vv*ny + nz*wv*wv*nz) / (3*r) ) / l_ref
                                                    + ( - (gam/(Re*Pr)) * E * r1 + (gam/(Re*Pr)) * 2.0*af * r1 ) / l_ref;
        
        tauDiff[i + (1+nd) * numPoints + 1 * sz2] = (1.0 / Re) * ( uv * r1 + nx*un / (3*r) ) / l_ref - ((gam/(Re*Pr)) * uv * r1) / l_ref;
        tauDiff[i + (1+nd) * numPoints + 2 * sz2] = (1.0 / Re) * ( vv * r1 + ny*un / (3*r) ) / l_ref - ((gam/(Re*Pr)) * vv * r1) / l_ref;
        tauDiff[i + (1+nd) * numPoints + 3 * sz2] = (1.0 / Re) * ( wv * r1 + nz*un / (3*r) ) / l_ref - ((gam/(Re*Pr)) * wv * r1) / l_ref;
        
        tauDiff[i + (1+nd) * numPoints + (1+nd) * sz2] = ((gam/(Re*Pr)) * r1) / l_ref;

        if (computeJacobian == 1) {
            r1m1 = -1.0 / (r * r);
            r1m2 = 0.0;
            r1m3 = 0.0;
            r1m4 = 0.0;
            r1m5 = 0.0;

            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;//r1
            uvm3 = ru * r1m3;//0.0
            uvm4 = ru * r1m4;//0.0
            uvm5 = ru * r1m5;//0.0

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;//0.0
            vvm3 = r1 + rv * r1m3;//r1
            vvm4 = rv * r1m4;//0.0
            vvm5 = rv * r1m5;//0.0
            
            wvm1 = rw * r1m1;
            wvm2 = rw * r1m2;//0.0
            wvm3 = rw * r1m3;//0.0
            wvm4 = r1 + rw * r1m4;//r1
            wvm5 = rw * r1m5;//0.0
            
            uvuvm1 = 2.0*uv*uvm1;
            uvuvm2 = 2.0*uv*uvm2;
            uvuvm3 = 2.0*uv*uvm3;//0.0
            uvuvm4 = 2.0*uv*uvm4;//0.0
            uvuvm5 = 2.0*uv*uvm5;//0.0
            
            uvvvm1 = uv*vvm1 + uvm1*vv;
            uvvvm2 = uv*vvm2 + uvm2*vv;
            uvvvm3 = uv*vvm3 + uvm3*vv;
            uvvvm4 = uv*vvm4 + uvm4*vv;//0.0
            uvvvm5 = uv*vvm5 + uvm5*vv;//0.0
            
            uvwvm1 = uv*wvm1 + uvm1*wv;
            uvwvm2 = uv*wvm2 + uvm2*wv;
            uvwvm3 = uv*wvm3 + uvm3*wv;//0.0
            uvwvm4 = uv*wvm4 + uvm4*wv;
            uvwvm5 = uv*wvm5 + uvm5*wv;//0.0
            
            vvuvm1 = vv*uvm1 + vvm1*uv;
            vvuvm2 = vv*uvm2 + vvm2*uv;
            vvuvm3 = vv*uvm3 + vvm3*uv;
            vvuvm4 = vv*uvm4 + vvm4*uv;//0.0
            vvuvm5 = vv*uvm5 + vvm5*uv;//0.0
            
            vvvvm1 = 2.0*vv*vvm1;
            vvvvm2 = 2.0*vv*vvm2;//0.0
            vvvvm3 = 2.0*vv*vvm3;
            vvvvm4 = 2.0*vv*vvm4;//0.0
            vvvvm5 = 2.0*vv*vvm5;//0.0
            
            vvwvm1 = vv*wvm1 + vvm1*wv;
            vvwvm2 = vv*wvm2 + vvm2*wv;//0.0
            vvwvm3 = vv*wvm3 + vvm3*wv;
            vvwvm4 = vv*wvm4 + vvm4*wv;
            vvwvm5 = vv*wvm5 + vvm5*wv;//0.0
            
            wvuvm1 = wv*uvm1 + wvm1*uv;
            wvuvm2 = wv*uvm2 + wvm2*uv;
            wvuvm3 = wv*uvm3 + wvm3*uv;//0.0
            wvuvm4 = wv*uvm4 + wvm4*uv;
            wvuvm5 = wv*uvm5 + wvm5*uv;//0.0
            
            wvvvm1 = wv*vvm1 + wvm1*vv;
            wvvvm2 = wv*vvm2 + wvm2*vv;//0.0
            wvvvm3 = wv*vvm3 + wvm3*vv;
            wvvvm4 = wv*vvm4 + wvm4*vv;
            wvvvm5 = wv*vvm5 + wvm5*vv;//0.0
            
            wvwvm1 = 2.0*wv*wvm1;
            wvwvm2 = 2.0*wv*wvm2;//0.0
            wvwvm3 = 2.0*wv*wvm3;//0.0
            wvwvm4 = 2.0*wv*wvm4;
            wvwvm5 = 2.0*wv*wvm5;//0.0
            
            uvuvr1m1 = uvuvm1*r1 + uvuv*r1m1;
            uvuvr1m2 = uvuvm2*r1 + uvuv*r1m2;
            uvuvr1m3 = uvuvm3*r1 + uvuv*r1m3;
            uvuvr1m4 = uvuvm4*r1 + uvuv*r1m4;
            uvuvr1m5 = uvuvm5*r1 + uvuv*r1m5;
            
            uvvvr1m1 = uvvvm1*r1 + uvvv*r1m1;
            uvvvr1m2 = uvvvm2*r1 + uvvv*r1m2;
            uvvvr1m3 = uvvvm3*r1 + uvvv*r1m3;
            uvvvr1m4 = uvvvm4*r1 + uvvv*r1m4;
            uvvvr1m5 = uvvvm5*r1 + uvvv*r1m5;
            
            uvwvr1m1 = uvwvm1*r1 + uvwv*r1m1;
            uvwvr1m2 = uvwvm2*r1 + uvwv*r1m2;
            uvwvr1m3 = uvwvm3*r1 + uvwv*r1m3;
            uvwvr1m4 = uvwvm4*r1 + uvwv*r1m4;
            uvwvr1m5 = uvwvm5*r1 + uvwv*r1m5;
            
            vvuvr1m1 = vvuvm1*r1 + vvuv*r1m1;
            vvuvr1m2 = vvuvm2*r1 + vvuv*r1m2;
            vvuvr1m3 = vvuvm3*r1 + vvuv*r1m3;
            vvuvr1m4 = vvuvm4*r1 + vvuv*r1m4;
            vvuvr1m5 = vvuvm5*r1 + vvuv*r1m5;
            
            vvvvr1m1 = vvvvm1*r1 + vvvv*r1m1;
            vvvvr1m2 = vvvvm2*r1 + vvvv*r1m2;
            vvvvr1m3 = vvvvm3*r1 + vvvv*r1m3;
            vvvvr1m4 = vvvvm4*r1 + vvvv*r1m4;
            vvvvr1m5 = vvvvm5*r1 + vvvv*r1m5;
            
            vvwvr1m1 = vvwvm1*r1 + vvwv*r1m1;
            vvwvr1m2 = vvwvm2*r1 + vvwv*r1m2;
            vvwvr1m3 = vvwvm3*r1 + vvwv*r1m3;
            vvwvr1m4 = vvwvm4*r1 + vvwv*r1m4;
            vvwvr1m5 = vvwvm5*r1 + vvwv*r1m5;
            
            wvuvr1m1 = wvuvm1*r1 + wvuv*r1m1;
            wvuvr1m2 = wvuvm2*r1 + wvuv*r1m2;
            wvuvr1m3 = wvuvm3*r1 + wvuv*r1m3;
            wvuvr1m4 = wvuvm4*r1 + wvuv*r1m4;
            wvuvr1m5 = wvuvm5*r1 + wvuv*r1m5;
            
            wvvvr1m1 = wvvvm1*r1 + wvvv*r1m1;
            wvvvr1m2 = wvvvm2*r1 + wvvv*r1m2;
            wvvvr1m3 = wvvvm3*r1 + wvvv*r1m3;
            wvvvr1m4 = wvvvm4*r1 + wvvv*r1m4;
            wvvvr1m5 = wvvvm5*r1 + wvvv*r1m5;
            
            wvwvr1m1 = wvwvm1*r1 + wvwv*r1m1;
            wvwvr1m2 = wvwvm2*r1 + wvwv*r1m2;
            wvwvr1m3 = wvwvm3*r1 + wvwv*r1m3;
            wvwvr1m4 = wvwvm4*r1 + wvwv*r1m4;
            wvwvr1m5 = wvwvm5*r1 + wvwv*r1m5;
            
            Em1 = rE * r1m1;
            Em2 = rE * r1m2;//0.0
            Em3 = rE * r1m3;//0.0
            Em4 = rE * r1m4;//0.0
            Em5 = r1 + rE * r1m5;

            afm1 = uv * uvm1 + vv * vvm1 + wv * wvm1;
            afm2 = uv * uvm2 + vv * vvm2 + wv * wvm2;
            afm3 = uv * uvm3 + vv * vvm3 + wv * wvm3;
            afm4 = uv * uvm4 + vv * vvm4 + wv * wvm4;
            afm5 = uv * uvm5 + vv * vvm5 + wv * wvm5;//0.0

            unm1 = uvm1 * nx + vvm1 * ny + wvm1 * nz;
            unm2 = uvm2 * nx + vvm2 * ny + wvm2 * nz;
            unm3 = uvm3 * nx + vvm3 * ny + wvm3 * nz;
            unm4 = uvm4 * nx + vvm4 * ny + wvm4 * nz;
            unm5 = uvm5 * nx + vvm5 * ny + wvm5 * nz;//0.0
            
            Er1m1 = Em1*r1+E*r1m1;
            Er1m2 = Em2*r1+E*r1m2;
            Er1m3 = Em3*r1+E*r1m3;
            Er1m4 = Em4*r1+E*r1m4;
            Er1m5 = Em5*r1+E*r1m5;
            
            uvr1m1 = uvm1*r1+uv*r1m1;
            uvr1m2 = uvm2*r1+uv*r1m2;
            uvr1m3 = uvm3*r1+uv*r1m3;//0.0
            uvr1m4 = uvm4*r1+uv*r1m4;//0.0
            uvr1m5 = uvm5*r1+uv*r1m5;//0.0
            
            vvr1m1 = vvm1*r1+vv*r1m1;
            vvr1m2 = vvm2*r1+vv*r1m2;//0.0
            vvr1m3 = vvm3*r1+vv*r1m3;
            vvr1m4 = vvm4*r1+vv*r1m4;//0.0
            vvr1m5 = vvm5*r1+vv*r1m5;//0.0
            
            wvr1m1 = wvm1*r1+wv*r1m1;
            wvr1m2 = wvm2*r1+wv*r1m2;//0.0
            wvr1m3 = wvm3*r1+wv*r1m3;//0.0
            wvr1m4 = wvm4*r1+wv*r1m4;
            wvr1m5 = wvm5*r1+wv*r1m5;//0.0
            
            unr1m1 = unm1*r1+un*r1m1;
            unr1m2 = unm2*r1+un*r1m2;
            unr1m3 = unm3*r1+un*r1m3;
            unr1m4 = unm4*r1+un*r1m4;
            unr1m5 = unm5*r1+un*r1m5;//0.0
            
            afr1m1 = afm1*r1+af*r1m1;
            afr1m2 = afm2*r1+af*r1m2;
            afr1m3 = afm3*r1+af*r1m3;
            afr1m4 = afm4*r1+af*r1m4;
            afr1m5 = afm5*r1+af*r1m5;//0.0
            
            tauDiff_uh[i + 1 * numPoints + 0 * sz2 + 0 * sz3] = - (1.0 / Re) * ( uvr1m1 + nx * unr1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 0 * sz2 + 1 * sz3] = - (1.0 / Re) * ( uvr1m2 + nx * unr1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 0 * sz2 + 2 * sz3] = - (1.0 / Re) * ( uvr1m3 + nx * unr1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 0 * sz2 + 3 * sz3] = - (1.0 / Re) * ( uvr1m4 + nx * unr1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 0 * sz2 + 4 * sz3] = - (1.0 / Re) * ( uvr1m5 + nx * unr1m5 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 2 * numPoints + 0 * sz2 + 0 * sz3] = - (1.0 / Re) * ( vvr1m1 + ny * unr1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 0 * sz2 + 1 * sz3] = - (1.0 / Re) * ( vvr1m2 + ny * unr1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 0 * sz2 + 2 * sz3] = - (1.0 / Re) * ( vvr1m3 + ny * unr1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 0 * sz2 + 3 * sz3] = - (1.0 / Re) * ( vvr1m4 + ny * unr1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 0 * sz2 + 4 * sz3] = - (1.0 / Re) * ( vvr1m5 + ny * unr1m5 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 3 * numPoints + 0 * sz2 + 0 * sz3] = - (1.0 / Re) * ( wvr1m1 + nz * unr1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 0 * sz2 + 1 * sz3] = - (1.0 / Re) * ( wvr1m2 + nz * unr1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 0 * sz2 + 2 * sz3] = - (1.0 / Re) * ( wvr1m3 + nz * unr1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 0 * sz2 + 3 * sz3] = - (1.0 / Re) * ( wvr1m4 + nz * unr1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 0 * sz2 + 4 * sz3] = - (1.0 / Re) * ( wvr1m5 + nz * unr1m5 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 1 * numPoints + 1 * sz2 + 0 * sz3] = (1.0 / Re) * ( r1m1 + nx * nx * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 1 * sz2 + 1 * sz3] = (1.0 / Re) * ( r1m2 + nx * nx * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 1 * sz2 + 2 * sz3] = (1.0 / Re) * ( r1m3 + nx * nx * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 1 * sz2 + 3 * sz3] = (1.0 / Re) * ( r1m4 + nx * nx * r1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 1 * sz2 + 4 * sz3] = (1.0 / Re) * ( r1m5 + nx * nx * r1m5 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 2 * numPoints + 1 * sz2 + 0 * sz3] = (1.0 / Re) * (           nx * ny * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 1 * sz2 + 1 * sz3] = (1.0 / Re) * (           nx * ny * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 1 * sz2 + 2 * sz3] = (1.0 / Re) * (           nx * ny * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 1 * sz2 + 3 * sz3] = (1.0 / Re) * (           nx * ny * r1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 1 * sz2 + 4 * sz3] = (1.0 / Re) * (           nx * ny * r1m5 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 3 * numPoints + 1 * sz2 + 0 * sz3] = (1.0 / Re) * (           nx * nz * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 1 * sz2 + 1 * sz3] = (1.0 / Re) * (           nx * nz * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 1 * sz2 + 2 * sz3] = (1.0 / Re) * (           nx * nz * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 1 * sz2 + 3 * sz3] = (1.0 / Re) * (           nx * nz * r1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 1 * sz2 + 4 * sz3] = (1.0 / Re) * (           nx * nz * r1m5 / 3.0 ) / l_ref;

            tauDiff_uh[i + 1 * numPoints + 2 * sz2 + 0 * sz3] = (1.0 / Re) * (           ny * nx * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 2 * sz2 + 1 * sz3] = (1.0 / Re) * (           ny * nx * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 2 * sz2 + 2 * sz3] = (1.0 / Re) * (           ny * nx * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 2 * sz2 + 3 * sz3] = (1.0 / Re) * (           ny * nx * r1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 2 * sz2 + 4 * sz3] = (1.0 / Re) * (           ny * nx * r1m5 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 2 * numPoints + 2 * sz2 + 0 * sz3] = (1.0 / Re) * ( r1m1 + ny * ny * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 2 * sz2 + 1 * sz3] = (1.0 / Re) * ( r1m2 + ny * ny * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 2 * sz2 + 2 * sz3] = (1.0 / Re) * ( r1m3 + ny * ny * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 2 * sz2 + 3 * sz3] = (1.0 / Re) * ( r1m4 + ny * ny * r1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 2 * sz2 + 4 * sz3] = (1.0 / Re) * ( r1m5 + ny * ny * r1m5 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 3 * numPoints + 2 * sz2 + 0 * sz3] = (1.0 / Re) * (           ny * nz * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 2 * sz2 + 1 * sz3] = (1.0 / Re) * (           ny * nz * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 2 * sz2 + 2 * sz3] = (1.0 / Re) * (           ny * nz * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 2 * sz2 + 3 * sz3] = (1.0 / Re) * (           ny * nz * r1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 2 * sz2 + 4 * sz3] = (1.0 / Re) * (           ny * nz * r1m5 / 3.0 ) / l_ref;

            tauDiff_uh[i + 1 * numPoints + 3 * sz2 + 0 * sz3] = (1.0 / Re) * (           nz * nx * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 3 * sz2 + 1 * sz3] = (1.0 / Re) * (           nz * nx * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 3 * sz2 + 2 * sz3] = (1.0 / Re) * (           nz * nx * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 3 * sz2 + 3 * sz3] = (1.0 / Re) * (           nz * nx * r1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 1 * numPoints + 3 * sz2 + 4 * sz3] = (1.0 / Re) * (           nz * nx * r1m5 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 2 * numPoints + 3 * sz2 + 0 * sz3] = (1.0 / Re) * (           nz * ny * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 3 * sz2 + 1 * sz3] = (1.0 / Re) * (           nz * ny * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 3 * sz2 + 2 * sz3] = (1.0 / Re) * (           nz * ny * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 3 * sz2 + 3 * sz3] = (1.0 / Re) * (           nz * ny * r1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 2 * numPoints + 3 * sz2 + 4 * sz3] = (1.0 / Re) * (           nz * ny * r1m5 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + 3 * numPoints + 3 * sz2 + 0 * sz3] = (1.0 / Re) * ( r1m1 + nz * nz * r1m1 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 3 * sz2 + 1 * sz3] = (1.0 / Re) * ( r1m2 + nz * nz * r1m2 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 3 * sz2 + 2 * sz3] = (1.0 / Re) * ( r1m3 + nz * nz * r1m3 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 3 * sz2 + 3 * sz3] = (1.0 / Re) * ( r1m4 + nz * nz * r1m4 / 3.0 ) / l_ref;
            tauDiff_uh[i + 3 * numPoints + 3 * sz2 + 4 * sz3] = (1.0 / Re) * ( r1m5 + nz * nz * r1m5 / 3.0 ) / l_ref;
            
            tauDiff_uh[i + (1+nd) * numPoints + 0 * sz2 + 0 * sz3] = - (1.0 / Re) * ( 2.0*afr1m1 + (nx*uvuvr1m1*nx + nx*uvvvr1m1*ny + nx*uvwvr1m1*nz + ny*vvuvr1m1*nx + ny*vvvvr1m1*ny + ny*vvwvr1m1*nz + nz*wvuvr1m1*nx + nz*wvvvr1m1*ny + nz*wvwvr1m1*nz) / 3.0 ) / l_ref
                                                                    + ( - (gam/(Re*Pr)) * Er1m1 + (gam/(Re*Pr)) * 2.0*afr1m1 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 0 * sz2 + 1 * sz3] = - (1.0 / Re) * ( 2.0*afr1m2 + (nx*uvuvr1m2*nx + nx*uvvvr1m2*ny + nx*uvwvr1m2*nz + ny*vvuvr1m2*nx + ny*vvvvr1m2*ny + ny*vvwvr1m2*nz + nz*wvuvr1m2*nx + nz*wvvvr1m2*ny + nz*wvwvr1m2*nz) / 3.0 ) / l_ref
                                                                    + ( - (gam/(Re*Pr)) * Er1m2 + (gam/(Re*Pr)) * 2.0*afr1m2 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 0 * sz2 + 2 * sz3] = - (1.0 / Re) * ( 2.0*afr1m3 + (nx*uvuvr1m3*nx + nx*uvvvr1m3*ny + nx*uvwvr1m3*nz + ny*vvuvr1m3*nx + ny*vvvvr1m3*ny + ny*vvwvr1m3*nz + nz*wvuvr1m3*nx + nz*wvvvr1m3*ny + nz*wvwvr1m3*nz) / 3.0 ) / l_ref
                                                                    + ( - (gam/(Re*Pr)) * Er1m3 + (gam/(Re*Pr)) * 2.0*afr1m3 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 0 * sz2 + 3 * sz3] = - (1.0 / Re) * ( 2.0*afr1m4 + (nx*uvuvr1m4*nx + nx*uvvvr1m4*ny + nx*uvwvr1m4*nz + ny*vvuvr1m4*nx + ny*vvvvr1m4*ny + ny*vvwvr1m4*nz + nz*wvuvr1m4*nx + nz*wvvvr1m4*ny + nz*wvwvr1m4*nz) / 3.0 ) / l_ref
                                                                    + ( - (gam/(Re*Pr)) * Er1m4 + (gam/(Re*Pr)) * 2.0*afr1m4 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 0 * sz2 + 4 * sz3] = - (1.0 / Re) * ( 2.0*afr1m5 + (nx*uvuvr1m5*nx + nx*uvvvr1m5*ny + nx*uvwvr1m5*nz + ny*vvuvr1m5*nx + ny*vvvvr1m5*ny + ny*vvwvr1m5*nz + nz*wvuvr1m5*nx + nz*wvvvr1m5*ny + nz*wvwvr1m5*nz) / 3.0 ) / l_ref
                                                                    + ( - (gam/(Re*Pr)) * Er1m5 + (gam/(Re*Pr)) * 2.0*afr1m5 ) / l_ref;
            
            tauDiff_uh[i + (1+nd) * numPoints + 1 * sz2 + 0 * sz3] = (1.0 / Re) * ( uvr1m1 + nx*unr1m1 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * uvr1m1 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 1 * sz2 + 1 * sz3] = (1.0 / Re) * ( uvr1m2 + nx*unr1m2 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * uvr1m2 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 1 * sz2 + 2 * sz3] = (1.0 / Re) * ( uvr1m3 + nx*unr1m3 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * uvr1m3 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 1 * sz2 + 3 * sz3] = (1.0 / Re) * ( uvr1m4 + nx*unr1m4 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * uvr1m4 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 1 * sz2 + 4 * sz3] = (1.0 / Re) * ( uvr1m5 + nx*unr1m5 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * uvr1m5 ) / l_ref;
            
            tauDiff_uh[i + (1+nd) * numPoints + 2 * sz2 + 0 * sz3] = (1.0 / Re) * ( vvr1m1 + ny*unr1m1 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * vvr1m1 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 2 * sz2 + 1 * sz3] = (1.0 / Re) * ( vvr1m2 + ny*unr1m2 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * vvr1m2 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 2 * sz2 + 2 * sz3] = (1.0 / Re) * ( vvr1m3 + ny*unr1m3 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * vvr1m3 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 2 * sz2 + 3 * sz3] = (1.0 / Re) * ( vvr1m4 + ny*unr1m4 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * vvr1m4 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 2 * sz2 + 4 * sz3] = (1.0 / Re) * ( vvr1m5 + ny*unr1m5 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * vvr1m5 ) / l_ref;
            
            tauDiff_uh[i + (1+nd) * numPoints + 3 * sz2 + 0 * sz3] = (1.0 / Re) * ( wvr1m1 + nz*unr1m1 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * wvr1m1 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 3 * sz2 + 1 * sz3] = (1.0 / Re) * ( wvr1m2 + nz*unr1m2 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * wvr1m2 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 3 * sz2 + 2 * sz3] = (1.0 / Re) * ( wvr1m3 + nz*unr1m3 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * wvr1m3 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 3 * sz2 + 3 * sz3] = (1.0 / Re) * ( wvr1m4 + nz*unr1m4 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * wvr1m4 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + 3 * sz2 + 4 * sz3] = (1.0 / Re) * ( wvr1m5 + nz*unr1m5 / 3.0 ) / l_ref - ((gam/(Re*Pr)) * wvr1m5 ) / l_ref;
            
            tauDiff_uh[i + (1+nd) * numPoints + (1+nd) * sz2 + 0 * sz3] = ((gam/(Re*Pr)) * r1m1 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + (1+nd) * sz2 + 1 * sz3] = ((gam/(Re*Pr)) * r1m2 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + (1+nd) * sz2 + 2 * sz3] = ((gam/(Re*Pr)) * r1m3 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + (1+nd) * sz2 + 3 * sz3] = ((gam/(Re*Pr)) * r1m4 ) / l_ref;
            tauDiff_uh[i + (1+nd) * numPoints + (1+nd) * sz2 + 4 * sz3] = ((gam/(Re*Pr)) * r1m5 ) / l_ref;
        }
    }
    // tauDiff: numPoints / nch / nch
    // tauDiff_uh: numPoints / nch / nch / nch
    
   for (i=0; i<sz4; i++)
       tauDiff_uh[i] = 0.0;
}

void getTensorialDiffStabilizationTensor(double* tauDiff, double* tauDiff_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, int numPoints, int nch, int nd, int computeJacobian)
{
    if (nd == 2)
        getTensorialDiffStabilizationTensor2d(tauDiff, tauDiff_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (nd == 3)
        getTensorialDiffStabilizationTensor3d(tauDiff, tauDiff_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else
        error("Number of dimensions not implemented.\n");
    // tau: numPoints / nch / nch
    // tau_uh: numPoints / nch / nch / nch
}

void getDiffusiveStabilizationTensor(double* tauDiff, double* tauDiff_uh, double* uh, double* pg, double* nl, appstruct &app, double* param, Int diffStabMethod, int numPoints, int nch, int nd, int computeJacobian)
{
    if (diffStabMethod == 0) {
    }
    if (diffStabMethod == 1)
        getConstantDiffStabilizationTensor(tauDiff, tauDiff_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else if (diffStabMethod == 2)
        getTensorialDiffStabilizationTensor(tauDiff, tauDiff_uh, uh, pg, nl, app, param, numPoints, nch, nd, computeJacobian);
    else
        error("Stabilization method for diffusive operator not implemented.\n");
}

void chainRuleJacobianStabilizationTensor(double* tau_UH, double* tau_uh, double* pg, appstruct &app, int numPoints, int nch, int nd)
{
    /* This function should be used only if ALEflag == 3 */

    int i, j, sz;

    Int ALEflag = app.ALEflag;

    /* Get gg */
    double * gg = &pg[(3 * nd + nd * nd) * numPoints];

    // Apply chain rule
    sz = nch * nch * nch;
    for (j = 0; j < sz; j++)
        for (i = 0; i < numPoints; i++)
            tau_UH[i + j * numPoints] = tau_uh[i + j * numPoints] / gg[i];
}

void getStabilizationTensor(double* tau, double* tau_UH, double* UH, double* pg, double* NL, appstruct &app, double* param,
                            int numPoints, int nch, int nd, int computeJacobian)
{
    Int inc = 1, len;
    int i;
    double one = 1.0;

    if (app.flag_q == 0)
        app.diffStabMethod = 0;

    Int convStabMethod = app.convStabMethod;
    Int diffStabMethod = app.diffStabMethod;
    Int ALEflag = app.ALEflag;

    // Get uh:
    double * uh;
    double * tau_uh;
    if (ALEflag == 3) {
        uh = new double[numPoints * nch];
        tau_uh = new double[numPoints * nch * nch * nch];
        UH2uh(uh, UH, pg, app, numPoints, nch, nd);
    }
    else {
        uh = &UH[0];
        tau_uh = &tau_UH[0];
    }

    // Get nl:
    double * nl;
    if (ALEflag == 2 || ALEflag == 3) {
        nl = new double[numPoints * nd];
        NL2nl(nl, NL, pg, app, numPoints, nd);
    }
    else {
        nl = &NL[0];
    }
    
    // Get stabilization tensor and derivatives in physical coordinates
    getConvectiveStabilizationTensor(tau, tau_uh, uh, pg, nl, app, param, convStabMethod, numPoints, nch, nd, computeJacobian);
    
    if (diffStabMethod != 0) {
        double *tauDiff  = new double [numPoints * nch * nch];
        double *tauDiff_uh  = new double [numPoints * nch * nch * nch];

        getDiffusiveStabilizationTensor(tauDiff, tauDiff_uh, uh, pg, nl, app, param, diffStabMethod, numPoints, nch, nd, computeJacobian);

        len = numPoints*nch*nch;
        DAXPY(&len, &one, &tauDiff[0], &inc, &tau[0], &inc);

        len = numPoints*nch*nch*nch;
        DAXPY(&len, &one, &tauDiff_uh[0], &inc, &tau_uh[0], &inc);

        delete[] tauDiff; delete[] tauDiff_uh;
    }

    // Apply chain rule
    if (ALEflag == 3 && computeJacobian == 1) {
        chainRuleJacobianStabilizationTensor(tau_UH, tau_uh, pg, app, numPoints, nch, nd);
    }
    
    // Deallocate dynamic memory
    if (ALEflag == 3) {
        delete[] uh; delete[] tau_uh;
    }
    if (ALEflag == 2 || ALEflag == 3) {
        delete[] nl;
    }
}

#endif
