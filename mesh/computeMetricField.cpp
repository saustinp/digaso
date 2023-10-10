#ifndef __COMPUTEMETRICFIELD
#define __COMPUTEMETRICFIELD

// Written by: C. Nguyen & P. Fernandez

void computeMetricField(meshstruct &mesh, masterstruct &master, Int *ndims, Int isotropicMasterFlag)
{
    Int i, j, k, l, ll, kA, kB;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int nve = ndims[3];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc = ndims[19];
    char chn = 'N', cht = 'T';
    double one = 1.0, zero = 0.0;
    
    Int elemtype = getElemtype(nd, nve);
    
    double *shapnvt = &master.shapnvt[0];       // npv (eval points) / npv (basis) / nd1
    double *pn;
    
    double *A = new double[nd*nd];
    double *J = new double[npv*nd*nd];      // npv / nd (x) / nd (xi)
    double *B = new double[nd*nd];
    double *M = new double[nd*nd];
    
    mesh.M.resize(ne*npv*nd*nd);
    
    if (isotropicMasterFlag == 0) {
        if (nd == 2) {
            A[0] = 1.0; A[2] = 0.0;
            A[1] = 0.0; A[3] = 1.0;
        }
        else if (nd == 3) {
            A[0] = 1.0; A[3] = 0.0; A[6] = 0.0;
            A[1] = 0.0; A[4] = 1.0; A[7] = 0.0;
            A[2] = 0.0; A[5] = 0.0; A[8] = 1.0;
        }
        else
            error("Number of dimensions not implemented.\n");
    }
    else if (isotropicMasterFlag == 1) {
        if (nd == 2 && elemtype == 0) {
            A[0] = 1.0; A[2] = -1.0/sqrt(3.0);
            A[1] = 0.0; A[3] = 2.0/sqrt(3.0);
        }
        else if (nd == 2 && elemtype == 1) {
            A[0] = 1.0; A[2] = 0.0;
            A[1] = 0.0; A[3] = 1.0;
        }  
        else if (nd == 3 && elemtype == 0) {
            A[0] = 1.0; A[3] = -sqrt(1.0/3.0); A[6] = -sqrt(1.0/6.0);
            A[1] = 0.0; A[4] =  sqrt(4.0/3.0); A[7] = -sqrt(1.0/6.0);
            A[2] = 0.0; A[5] =  0.0;           A[8] =  sqrt(3.0/2.0);
        }
        else if (nd == 3 && elemtype == 1) {
            A[0] = 1.0; A[3] = 0.0; A[6] = 0.0;
            A[1] = 0.0; A[4] = 1.0; A[7] = 0.0;
            A[2] = 0.0; A[5] = 0.0; A[8] = 1.0;
        }
        else
            error("Element type or number of dimensions not implemented.\n");
    }
    
    for (i = 0; i < ne; i++) {
        pn = &mesh.dgnodes[i*ncd*npv];      // npv / nd
        
        // Compute Jacobian of mapping from (isotropic?) master element to physical domain
        for (j=0; j<npv; j++) {
            for (k=0; k<nd; k++) {
                kA = (k+1)*npv*npv;
                kB = k*npv*nd;
                DGEMM(&chn, &chn, &npv, &nd, &npv, &one, &shapnvt[kA], &npv, &pn[0],
                      &npv, &zero, &J[kB], &npv);
            }
        }
        
        for (j = 0; j < npv; j++) {
            // Compute Jacobian of mapping from isotropic mater element to physical domain: B = J*A
            for (k = 0; k < nd; k++)
                for (l = 0; l < nd; l++) {
                    B[k*nd+l] = 0.0;
                    for (ll = 0; ll < nd; ll++)
                        B[k*nd+l] += J[ll*npv*nd+l*npv+j] * A[k*nd+ll];
                }
            
            // Compute inverse of metric tensor: M = B'*B
            DGEMM(&cht, &chn, &nd, &nd, &nd, &one, &B[0], &nd, &B[0],
                  &nd, &zero, &M[0], &nd);
            
            for (k = 0; k < nd*nd; k++)
                mesh.M[i*npv*nd*nd+j*nd*nd+k] = M[k];
        }
    }
    
    delete[] A; delete[] J;
    delete[] B; delete[] M;
}

void computeInvMetricField(meshstruct &mesh, masterstruct &master, Int *ndims, Int isotropicMasterFlag)
{
    Int i, j, k, l, ll, kA, kB, info;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int nve = ndims[3];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc = ndims[19];
    Int lwork = nd;
    char chn = 'N', cht = 'T';
    double one = 1.0, zero = 0.0;
    
    Int elemtype = getElemtype(nd, nve);
    
    double *shapnvt = &master.shapnvt[0];       // npv (eval points) / npv (basis) / nd1
    double *pn;
    
    Int *ipiv = new Int[nd];
    double *work = new double[nd];
    
    double *Binv; 
    double *A = new double[nd*nd];
    double *J = new double[npv*nd*nd];      // npv / nd (x) / nd (xi)
    double *B = new double[nd*nd];
    double *Minv = new double[nd*nd];
    
    mesh.Minv.resize(ne*npv*nd*nd);
    
    if (isotropicMasterFlag == 0) {
        if (nd == 2) {
            A[0] = 1.0; A[2] = 0.0;
            A[1] = 0.0; A[3] = 1.0;
        }
        else if (nd == 3) {
            A[0] = 1.0; A[3] = 0.0; A[6] = 0.0;
            A[1] = 0.0; A[4] = 1.0; A[7] = 0.0;
            A[2] = 0.0; A[5] = 0.0; A[8] = 1.0;
        }
        else
            error("Number of dimensions not implemented.\n");
    }
    else if (isotropicMasterFlag == 1) {
        if (nd == 2 && elemtype == 0) {
            A[0] = 1.0; A[2] = -1.0/sqrt(3.0);
            A[1] = 0.0; A[3] = 2.0/sqrt(3.0);
        }
        else if (nd == 2 && elemtype == 1) {
            A[0] = 1.0; A[2] = 0.0;
            A[1] = 0.0; A[3] = 1.0;
        }  
        else if (nd == 3 && elemtype == 0) {
            A[0] = 1.0; A[3] = -sqrt(1.0/3.0); A[6] = -sqrt(1.0/6.0);
            A[1] = 0.0; A[4] =  sqrt(4.0/3.0); A[7] = -sqrt(1.0/6.0);
            A[2] = 0.0; A[5] =  0.0;           A[8] =  sqrt(3.0/2.0);
        }
        else if (nd == 3 && elemtype == 1) {
            A[0] = 1.0; A[3] = 0.0; A[6] = 0.0;
            A[1] = 0.0; A[4] = 1.0; A[7] = 0.0;
            A[2] = 0.0; A[5] = 0.0; A[8] = 1.0;
        }
        else
            error("Element type or number of dimensions not implemented.\n");
    }
    
    for (i = 0; i < ne; i++) {
        pn = &mesh.dgnodes[i*ncd*npv];      // npv / nd
        
        // Compute Jacobian of mapping from (isotropic?) master element to physical domain
        for (j=0; j<npv; j++) {
            for (k=0; k<nd; k++) {
                kA = (k+1)*npv*npv;
                kB = k*npv*nd;
                DGEMM(&chn, &chn, &npv, &nd, &npv, &one, &shapnvt[kA], &npv, &pn[0],
                      &npv, &zero, &J[kB], &npv);
            }
        }
        
        for (j=0; j<npv; j++) {
            // Compute Jacobian of mapping from isotropic mater element to physical domain: B = J*A
            for (k = 0; k < nd; k++)
                for (l = 0; l < nd; l++) {
                    B[k*nd+l] = 0.0;
                    for (ll = 0; ll < nd; ll++)
                        B[k*nd+l] += J[ll*npv*nd+l*npv+j] * A[k*nd+ll];
                }
            
            // Compute inv(B):
            DGETRF(&nd, &nd, &B[0], &nd, &ipiv[0], &info);
            Binv = &B[0];
            DGETRI(&nd, &Binv[0], &nd, &ipiv[0], &work[0], &lwork, &info);
            
            // Compute inverse of metric tensor: Minv = Binv'*Binv
            DGEMM(&cht, &chn, &nd, &nd, &nd, &one, &Binv[0], &nd, &Binv[0],
                  &nd, &zero, &Minv[0], &nd);

            for (k = 0; k < nd*nd; k++)
                mesh.Minv[i*npv*nd*nd+j*nd*nd+k] = Minv[k];
        }
    }
    
    delete[] ipiv; delete[] work;
    delete[] A; delete[] J;
    delete[] B; delete[] Minv;
}

#endif
