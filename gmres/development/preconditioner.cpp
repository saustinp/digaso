#ifndef __PRECONDITIONER
#define __PRECONDITIONER

#include "ordering.cpp"

// Written by: C. Nguyen & P. Fernandez

Int sys.subdomainPreconditioner; // Type of subdomain preconditioner. -1: No preconditioner. 0: BILU(k). 1: Entity-based block Jacobi
Int sys.parallelPreconditioner; // Type of parallel preconditioner. Flags are TBD.
Int sys.numRowsP: // Number of entity-rows in subdomain preconditioner
Int sys.numRowsJ: // Number of entity-rows in subdomain Jacobian
Int sys.numRowsR; // Number of entity-rows in subdomain residual
Int sys.numBlocksP; // Number of blocks in subdomain preconditioner
Int sys.numBlocksJ; // Number of blocks in subdomain Jacobian
Int sys.numIntEnt;
vector <Int> sys.entSz[numEntities];
vector <Int> sys.ent2entP[numBlocksP];
vector <Int> sys.ent2entJ[numBlocksJ];
vector <Int> sys.ent2entStartP[numRowsP+1];
vector <Int> sys.ent2entStartJ[numRowsJ+1];
vector <Int> sys.rowStartP[numRowsP+1];
vector <Int> sys.rowStartJ[numRowsJ+1];
vector <Int> sys.rowStartR[numRowsR+1];
vector <Int> sys.blkStartP[numBlocksP+1];
vector <Int> sys.blkStartJ[numBlocksJ+1];
vector <Int> sys.MgMinusHg[numBlocksP]; // Entries of the spartsity pattern of Mg  that are also in Hg. Only for BILU(k) subdomain preconditioner.
        
void postproJPR(sysstruct &sys) {
    
    // Outputs:
    // - sys.rowStartJ, sys.rowStartP, sys.rowStartR
    // - sys.blkStartJ, sys.blkStartP
    
    // Inputs:
    // - sys.ent2entStartJ, sys.ent2entStartP
    // - sys.ent2entJ, sys.ent2entP
    // - sys.entSz
    
    Int numRowsJ = sys.numRowsJ;
    Int numRowsP = sys.numRowsP;
    Int numRowsR = sys.numRowsR;
    Int numBlocksJ = sys.numBlocksJ;
    Int numBlocksP = sys.numBlocksP;
    Int *rowStartJ = &sys.rowStartJ[0];
    Int *rowStartP = &sys.rowStartP[0];
    Int *rowStartR = &sys.rowStartR[0];
    Int *blkStartJ = &sys.blkStartJ[0];
    Int *blkStartP = &sys.blkStartP[0];
    Int *ent2entStartJ = &sys.ent2entStartJ[0];
    Int *ent2entStartP = &sys.ent2entStartP[0];
    Int *ent2entJ = &sys.ent2entJ[0];
    Int *ent2entP = &sys.ent2entP[0];
    Int *entSz = &sys.entSz[0];
    
    // Jacobian and residual:
    sys.rowStartJ.resize(numRowsJ);
    sys.rowStartJ[0] = 0;
    sys.rowStartR.resize(numRowsJ);
    sys.rowStartR[0] = 0;
    sys.blkStartJ.resize(numBlocksJ);
    sys.blkStartJ[0] = 0;
    Int i_b = 0;
    for (Int ri = 0; ri < numRowsJ; ri++) {
        Int piSp = sys.ent2entStartJ[ri];
        Int nri = sys.ent2entStartJ[ri+1] - sys.ent2entStartJ[ri];
        Int bsz_row = entSz[ri];
        Int bsz2s = 0;
        
        for (Int j = 0; j < nri; j++) {
            Int rj = ent2ent[piSp+j];
            Int bsz_col = entSz[rj];
            Int bsz2 = bsz_row * bsz_col;
            sys.blkStartJ[i_b+1] = sys.blkStartJ[i_b] + bsz2;
            i_b++;
            bsz2s += bsz2;
        }
        sys.rowStartJ[ri+1] = sys.rowStartJ[ri] + bsz2s;
        sys.rowStartR[ri+1] = sys.rowStartR[ri] + bsz_row;
    }
    
    // Preconditioner:
    sys.rowStartP.resize(numRowsP);
    sys.rowStartP[0] = 0;
    sys.blkStartP.resize(numBlocksP);
    sys.blkStartP[0] = 0;
    i_b = 0;
    for (Int ri = 0; ri < numRowsP; ri++) {
        Int piSp = sys.ent2entStartP[ri];
        Int nri = sys.ent2entStartP[ri+1] - sys.ent2entStartP[ri];
        Int bsz_row = entSz[ri];
        Int bsz2s = 0;
        
        for (Int j = 0; j < nri; j++) {
            Int rj = ent2ent[piSp+j];
            Int bsz_col = entSz[rj];
            Int bsz2 = bsz_row * bsz_col;
            sys.blkStartP[i_b+1] = sys.blkStartP[i_b] + bsz2;
            i_b++;
            bsz2s += bsz2;
        }
        sys.rowStartP[ri+1] = sys.rowStartP[ri] + bsz2s;
    }
}

Int matchSparsity(sysstruct &sys, Int ri, Int j)
{
    Int piP = sys.ent2entStartP[ri];
    Int piJ = sys.ent2entStartJ[ri];
    Int *ent2entStartJ = &sys.ent2entStartJ[0];
    Int rj = sys.ent2entP[piP+j];
    Int *ent2entJ = &sys.ent2entJ[0];
    
    Int match = 0;
    Int nriJ = ent2entStartJ[ri+1] - ent2entStartJ[ri];
    for (Int k = 0; k < nriJ; k++) {
        if (ent2entJ[piJ+k] == rj) {
            match = 1;
            break;
        }
    }
    
    return match;
}

void computeMgMinusHg(sys)
{
    // Compute what entries of the spartsity pattern of Mg  
    // are also in Hg. The solution in stored in sys.MgMinusHg.
    // Only valid for BILU(k) subdomain preconditioner.
    
    if (sys.subdomainPreconditioner != 0)
        error("computeMgMinusHg only valid for BILU(k) subdomain preconditioner.\n");
    
    Int *ent2entStartP = &sys.ent2entStartP[0];
    
    sys.MgMinusHg.resize(sys.ent2entP.size());
    
    Int numRowsP = sys.numRowsP;
    Int count = 0;
    for (Int ri = 0; ri < numRowsP; ri++) {
        Int nriP = ent2entStartP[ri+1] - ent2entStartP[ri];
        for (Int i = 0; i < nriP; i++)
            sys.MgMinusHg[count++] = matchSparsity(sys, ri, i);
    }
    
    if (count != sys.ent2entP.size())
        error("Error H4VCB0 in computeMgMinusHg.\n");
}

void copyHg2MgWithMgFormat(sysstruct &sys)
{
    Int inc = 1;
    Int nrowsP = sys.numRowsP;
    
    Int *blkStartJ = &sys.blkStartJ[0];
    Int *blkStartP = &sys.blkStartP[0];
    Int *rowStartJ = &sys.rowStartJ[0];
    Int *rowStartP = &sys.rowStartP[0];
    Int *ent2entStartJ = &sys.ent2entStartJ[0];
    Int *ent2entStartP = &sys.ent2entStartP[0];
    Int *ent2entP = &sys.ent2entP[0];
    Int *Loent2entStart = &sys.Loent2entStart[0];
    Int *Uoent2entStart = &sys.Uoent2entStart[0];
    
    for (Int rjj=0; rjj<nrowsP; rjj++) {
        Int rj = sys.unordered2ordered[rjj];
        
        Int bsz_row = sys.entSz[rjj];
        Int bsz2 = bsz_row*bsz_row;
        
        Int pjP   = sys.UorowStart[rj]-bsz2;
        Int pjP_l = sys.LorowStart[rj];
        Int pjP_u = sys.UorowStart[rj+1];
        
        Int pjjJ  = rowStartJ[rjj];
        Int pjjP  = rowStartP[rjj];
        Int pjjSpP = ent2entStartP[rjj];
        
        Int nrjJ  = ent2entStartJ[rjj+1]-ent2entStartJ[rjj];
        Int nrjP  = ent2entStartP[rjj+1]-ent2entStartP[rjj];
        Int nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];        // Neighbors in L
        Int nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];        // Neighbors in U (includes the diagonal block)
        
        // Copy diagonal block:
        DCOPY(&bsz2, &sys.Hg[pjjJ], &inc, &sys.Mg[pjP], &inc);
        
        // Copy off-diagonal blocks:
        Int i_u = 0, i_l = 0;
        Int pjiP_l = pjP_l;
        Int pjiP_u = pjP_u;
        Int pjjiJ = pjjJ + bsz2;
        for (Int i=1; i<nrjP; i++) {
            rii = ent2entP[pjjSpP+i];
            ri = sys.unordered2ordered[rii];
            
            Int bsz_col = sys.entSz[rii];
            bsz2 = bsz_row*bsz_col;
            
            if (ri > rj) {
                // U block
                i_u ++;
                if (! sys.MgMinusHg[ent2entStartP[rjj]+i])
                    std::fill(&sys.Mg[pjiP_u], &sys.Mg[pjiP_u+bsz2], 0.0);
                else {
                    DCOPY(&bsz2, &sys.Hg[pjjiJ], &inc, &sys.Mg[pjiP_u], &inc);
                    pjjiJ += bsz2;
                }
                pji_u += bsz2;
            }
            else if (ri < rj) {
                // L block
                i_l ++;
                if (! sys.MgMinusHg[ent2entStartP[rjj]+i])
                    std::fill(&sys.Mg[pjiP_l], &sys.Mg[pjiP_l+bsz2], 0.0);
                else {
                    DCOPY(&bsz2, &sys.Hg[pjjiJ], &inc, &sys.Mg[pjiP_l], &inc);
                    pjjiJ += bsz2;
                }
                pji_l += bsz2;
            }
            else
                error("Error No. 1 in copyHg2MgWithMgFormat.\n");
        }
        
        if ((i_l != nrj_l || i_u != (nrj_u-1)) && nrjP > 0) {
            printf("Error No. 2 in copyHg2MgWithMgFormat. rjj= %d, rj = %d, i_l = %d, nrj_l = %d, i_u = %d, nrj_u = %d, my_rank = %d\n.", rjj, rj, i_l, nrj_l, i_u, nrj_u, sys.my_rank);
            error("\n");
        }
    }
}

void computeBILUk(sysstruct &sys)
{
    Int inc = 1, info;
    Int nrowsP = sys.numRowsP;
    double one = 1.0, zero = 0.0, minusone = -1.0;
    char chn = 'N';
    
    Int *ent2entStartJ = &sys.ent2entStartJ[0];
    Int *ent2entStartP = &sys.ent2entStartP[0];
    Int *ent2entJ = &sys.ent2entJ[0];
    Int *ent2entP = &sys.ent2entP[0];
    Int *blkStartJ = &sys.blkStartJ[0];
    Int *blkStartP = &sys.blkStartP[0];
    Int *rowStartJ = &sys.rowStartJ[0];
    Int *rowStartP = &sys.rowStartP[0];
    Int *entSz = &sys.entSz[0];
    Int *ipiv = &sys.ipiv[0];
    double *work = &sys.work[0];
    double *Lij = &sys.HiiInv[0];
    
    // Copy Hg into Mg with the desired storage format
    copyHg2MgWithMgFormat(sys);
    
    // Compute BILUk preconditioner
    for (Int rj=0; rj<nrowsP; rj++) {
        Int rjj = sys.ordered2unordered[rj];
        
        Int bsz_jj = entSz[rjj];
        Int bsz2 = bsz_jj*bsz_jj;
        Int lwork = bsz_jj;
        
        Int pjSp   = sys.Uoent2entStart[rj]-1;
        Int pjSp_l = sys.Loent2entStart[rj];
        Int pjSp_u = sys.Uoent2entStart[rj+1];
        Int pj   = sys.UorowStart[rj]-bsz2;
        Int pj_l = sys.LorowStart[rj];
        Int pj_u = sys.UorowStart[rj+1];
        
        // Compute inverse of block diagonal of Hg using LU decomposition
        DGETRF(&bsz_jj,&bsz_jj,&sys.Mg[pj],&bsz_jj,ipiv,&info);
        DGETRI(&bsz_jj,&sys.Mg[pj],&bsz_jj,ipiv,work,&lwork,&info);
        
        Int nrjP = ent2entStartP[rjj+1]-ent2entStartP[rjj];
        Int nrj_l = sys.Loent2entStart[rj+1]-sys.Loent2entStart[rj];        // Neighbors in L
        Int nrj_u = sys.Uoent2entStart[rj]-sys.Uoent2entStart[rj+1];        // Neighbors in U (including the diagonal block)
        for (Int i=0; i<(nrj_u-1); i++) {
            Int pjiSp_u = pjSp_u+i;
            Int pji_u = sys.LUoblkStart[pjiSp_u];
            Int rii = sys.LUoent2ent[pjiSp_u];
            Int ri = sys.unordered2ordered[rii];
            
            Int bsz_ii = sys.entSz[rii];
            bsz2 = bsz_jj*bsz_ii;
            
            Int piSp   = sys.Uoent2entStart[ri]-1;
            Int piSp_l = sys.Loent2entStart[ri];
            Int piSp_u = sys.Uoent2entStart[ri+1];
            Int pi   = sys.UorowStart[ri]-bsz_ii*bsz_ii;
            Int pi_l = sys.LorowStart[ri];
            Int pi_u = sys.UorowStart[ri+1]; 
                    
            Int nriP = ent2entStartP[rii+1]-ent2entStartP[rii];
            Int nri_l = sys.Loent2entStart[ri+1]-sys.Loent2entStart[ri];
            Int nri_u = sys.Uoent2entStart[ri]-sys.Uoent2entStart[ri+1];
            Int k;
            for (k=0; k<nri_l; k++) {
                if (sys.LUoent2ent[piSp_l+k] == rjj) break;
            }
            if (k==nri_l && k!=0) {
                error("Error No. 1 in computeBILUk.\n");
            }
            Int pijSp_l = piSp_l+k;
            Int pij_l = sys.LUoblkStart[pijSp_l];
            
            // Lij <- Uij*inv(Ujj)
            DGEMM(&chn, &chn, &bsz_ii, &bsz_jj, &bsz_jj, &one, &sys.Mg[pij_l], &bsz_ii, &sys.Mg[pj],
                  &bsz_jj, &zero, Lij, &bsz_ii);
            DCOPY(&bsz2, Lij, &inc, &sys.Mg[pij_l], &inc);
            
            // Uii <- Uii - Lij*Uji
            DGEMM(&chn, &chn, &bsz_ii, &bsz_jj, &bsz_ii, &minusone, &sys.Mg[pij_l], &bsz_ii, &sys.Mg[pji_u],
                  &bsz_jj, &one, &sys.Mg[pi], &bsz_ii);
            
            // L_il <- L_il - L_ij*U_jm (m "=" l)
            for (Int l=0; l<nri_l; l++) {
                Int pilSp_l = piSp_l+l;
                Int pil_l = sys.LUoblkStart[pilSp_l];
                Int rll  = sys.LUoent2ent[pilSp_l];
                Int rl = sys.unordered2ordered[rll];
                Int bsz_ll = entSz[rll];
                
                if (rl > rj) {
                    for (Int m=0; m<(nrj_u-1); m++) {
                        if (rll == sys.LUoent2ent[pjSp_u+m]) {
                            Int pjm_u = sys.LUoblkStart[pjSp_u+m];
                            DGEMM(&chn, &chn, &bsz_ii, &bsz_ll, &bsz_jj, &minusone, &sys.Mg[pij_l], &bsz_ii, &sys.Mg[pjm_u],
                                  &bsz_jj, &one, &sys.Mg[pil_l], &bsz_ii);
                        }
                    }
                }
            }

            // U_il <- U_il - L_ij*U_jm (m "=" l)
            for (Int l=0; l<(nri_u-1); l++) {
                Int pilSp_u = piSp_u+l;
                Int pil_u = sys.LUoblkStart[pilSp_u];
                Int rll  = sys.LUoent2ent[pilSp_u];
                Int rl = sys.unordered2ordered[rll];
                Int bsz_ll = entSz[rll];

                if (rl > rj) {
                    for (Int m=0; m<(nrj_u-1); m++) {
                        if (rll == sys.LUoent2ent[pjSp_u+m]) {
                            Int pjm_u = sys.LUoblkStart[pjSp_u+m];
                            DGEMM(&chn, &chn, &bsz_ii, &bsz_ll, &bsz_jj, &minusone, &sys.Mg[pij_l], &bsz_ii, &sys.Mg[pjm_u],
                                  &bsz_jj, &one, &sys.Mg[pil_u], &bsz_ii);
                        }
                    }
                }
            }
        }
    }
}

void computeEntityInvBJ(sysstruct &sys)
{
    Int info, inc = 1;
    double one = 1.0, zero = 0.0, minusone = -1.0;
    char chn = 'N';
    
    Int nrowsP = sys.numRowsP;
    Int *rowStartJ = &sys.rowStartJ[0];
    Int *rowStartP = &sys.rowStartP[0];
    
    vector<vector<Int> >* ipivs = &sys.ipivs;
    vector<vector<double> >* works = &sys.works;
    
//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for private(pj, info)
        for (Int rj=0; rj<nrowsP; rj++) {
            Int pjJ = rowStartJ[rj];
            Int pjP = rowStartP[rj];
            
            Int bsz = sys.entSz[rj];
            Int bsz2 = bsz*bsz;
            Int lwork = bsz;
            
            // Copy diagonal block:
            DCOPY(&bsz2, &sys.Hg[pjJ], &inc, &sys.Mg[pjP], &inc);
            
            // Compute inverse of block diagonal of Hg using LU decomposition:
            DGETRF(&bsz,&bsz,&sys.Mg[pjP],&bsz,&ipivs[0][this_thread][0],&info);
            DGETRI(&bsz,&sys.Mg[pjP],&bsz,&ipivs[0][this_thread][0],&works[0][this_thread][0],&lwork,&info);
        }
//     }
}

void applyBILUk_sp(sysstruct &sys, float * r_sp, double * r, Int * requestCounter)
{
    // Mg is stored in memory in the same order as is required for the BILUk solve (L first, then U).
    
    error("To be coded.\n");
    
#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && sys.parallelPreconditioner == ?? && (sys.reorderMethod == 0 || sys.reorderMethod == 1 || sys.reorderMethod == 2))
        recvvector(sys, &r[0], requestCounter);
#endif
    
    if (sys.reorderMethod == 3 || sys.reorderMethod == 4)
        error("Mixed-precision algorithm not implemented for MDF with constraints.\n");
    
    Int inc = 1, info;
    Int nrowsP = sys.numRowsP;
    Int numIntRows = sys.numIntRows;
    float one = 1.0, zero = 0.0, minusone = -1.0;
    char chn = 'N';
    
    float *rDenseRow_sp = &sys.rDenseRow_sp[0];
    Int *ordered2unordered = &sys.ordered2unordered[0];
    Int *Uoent2entStart = &sys.Uoent2entStart[0];
    Int *Loent2entStart = &sys.Loent2entStart[0];
    Int *LUoent2ent = &sys.LUoent2ent[0];
    Int *UorowStart = &sys.UorowStart[0];
    Int *LorowStart = &sys.LorowStart[0];
    Int *rowStartR = &sys.rowStartR[0];
    Int *entSz = &sy.entSz[0];
    
//     for (i = 0; i < nrowsP*bsz; i++)
//         r_sp[i] = (float) r[i];
    sys.r_sp.insert(sys.r_sp.begin(), sys.r.begin(), sys.r.end());
    
    // Forward solve r <- L \ r
    if (sys.precSolveImplementation == 0) {
        for (Int rj=0; rj<numIntRows; rj++) {
            Int rjj = ordered2unordered[rj];
            Int pj_l = LorowStart[rj];
            Int pjSp_l = Loent2entStart[rj];
            Int nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];
            
            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            Int pji_l = pj_l;
            for (Int i=0; i<nrj_l; i++) {
                Int rii = LUoent2ent[pjSp_l+i];
                Int piiR = rowStartR[rii];
                Int bsz_col = entSz[rii];
                Int bsz2 = bsz_row*bsz_col;
                
                SGEMV(&chn, &bsz_row, &bsz_col, &minusone, &sys.Mg_sp[pji_l], &bsz_row, &r_sp[piiR],
                      &inc, &one, &r_sp[pjjR], &inc);
                pji_l += bsz2;
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (Int rj=0; rj<numIntRows; rj++) {
            Int rjj = ordered2unordered[rj];
            Int pj_l = LorowStart[rj];
            Int pjSp_l = Loent2entStart[rj];
            Int nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            Int bsz_cols = 0;
            for (Int i=0; i<nrj_l; i++) {
                Int rii = LUoent2ent[pjSp_l+i];
                Int piiR = rowStartR[rii];
                Int bsz_col = entSz[rii];
                
                for (Int jjj = 0; jjj < bsz_col; jjj++)
                   rDenseRow_sp[bsz_cols+jjj] = r_sp[piiR+jjj];
                //SCOPY(&bsz, &r_sp[piiR], &inc, &rDenseRow_sp[bsz_cols], &inc);
                bsz_cols += bsz_col;
            }
            
            if (nrj_l > 0) {
                SGEMV(&chn, &bsz_row, &bsz_cols, &minusone, &sys.Mg_sp[pj_l], &bsz_row, &rDenseRow_sp[0],
                      &inc, &one, &r_sp[pjjR], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    
#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && parallelPreconditioner == ?? && (sys.reorderMethod == 3 || sys.reorderMethod == 4))
        recvvector(sys, &r[0], requestCounter);
#endif
    
    if (sys.precSolveImplementation == 0) {
        for (Int rj=numIntRows; rj<nrowsP; rj++) {
            Int rjj = ordered2unordered[rj];
            Int pj_l = LorowStart[rj];
            Int pjSp_l = Loent2entStart[rj];
            Int nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];
            
            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            Int pji_l = pj_l;
            for (Int i=0; i<nrj_l; i++) {
                Int rii = LUoent2ent[pjSp_l+i];
                Int piiR = rowStartR[rii];
                Int bsz_col = entSz[rii];
                Int bsz2 = bsz_row*bsz_col;
                
                SGEMV(&chn, &bsz_row, &bsz_col, &minusone, &sys.Mg_sp[pji_l], &bsz_row, &r_sp[piiR],
                      &inc, &one, &r_sp[pjjR], &inc);
                pji_l += bsz2;
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (Int rj=numIntRows; rj<nrowsP; rj++) {
            Int rjj = ordered2unordered[rj];
            Int pj_l = LorowStart[rj];
            Int pjSp_l = Loent2entStart[rj];
            Int nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            Int bsz_cols = 0;
            for (Int i=0; i<nrj_l; i++) {
                Int rii = LUoent2ent[pjSp_l+i];
                Int piiR = rowStartR[rii];
                Int bsz_col = entSz[rii];
                
                for (Int jjj = 0; jjj < bsz_col; jjj++)
                   rDenseRow_sp[bsz_cols+jjj] = r_sp[piiR+jjj];
                //SCOPY(&bsz, &r_sp[piiR], &inc, &rDenseRow_sp[bsz_cols], &inc);
                bsz_cols += bsz_col;
            }
            
            if (nrj_l > 0) {
                SGEMV(&chn, &bsz_row, &bsz_cols, &minusone, &sys.Mg_sp[pj_l], &bsz_row, &rDenseRow_sp[0],
                      &inc, &one, &r_sp[pjjR], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    
    // Backward solve r <- U \ r
    if (sys.precSolveImplementation == 0) {
        for (Int j=0; j<nrowsP; j++) {
            Int rj = nrowsP-j-1;
            Int rjj = ordered2unordered[rj];
            //Int pjSp = Uoent2entStart[rj]-1;
            Int pjSp_u = Uoent2entStart[rj+1];
            Int pj_u = UorowStart[rj+1];
            Int nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];
            
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];
            
            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            Int pji_u = pj_u;
            for (Int i=0; i<(nrj_u-1); i++) {
                Int rii = LUoent2ent[pjSp_u+i];
                Int piiR = rowStartR[rii];
                Int bsz_col = entSz[rii];
                
                SGEMV(&chn, &bsz_row, &bsz_col, &minusone, &sys.Mg_sp[pji_u], &bsz_row, &r_sp[piiR],
                      &inc, &one, &r_sp[pjjR], &inc);
                pji_u += bsz_row * bsz_col;
            }
            
            // Compute r_j <- U_jj \ r_j
            SGEMV(&chn, &bsz_row, &bsz_row, &one, &sys.Mg_sp[pji_u], &bsz_row, &r_sp[pjjR],
                  &inc, &zero, &sys.rtmp_sp[0], &inc);
            
            SCOPY(&bsz_row, &sys.rtmp_sp[0], &inc, &r_sp[pjjR], &inc);
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (Int j=0; j<nrowsP; j++) {
            Int rj = nrowsP-j-1;
            Int rjj = ordered2unordered[rj];
            //Int pjSp = Uoent2entStart[rj]-1;
            Int pjSp_u = Uoent2entStart[rj+1];
            Int pj_u = UorowStart[rj+1];
            Int nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];
            
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];
            
            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            Int pji_u = pj_u;
            Int bsz_cols = 0;
            for (Int i=0; i<(nrj_u-1); i++) {
                Int rii = LUoent2ent[pjSp_u+i];
                Int piiR = rowStartR[rii];
                Int bsz_col = entSz[rii];
                Int bsz2 = bsz_row * bsz_col;
                Int pjiSp_u = pjSp_u+i;
                
                for (Int jjj = 0; jjj < bsz_col; jjj++)
                    rDenseRow_sp[bsz_cols+jjj] = r_sp[piiR+jjj];
                //SCOPY(&bsz_col, &r_sp[piiR], &inc, &rDenseRow_sp[bsz_cols], &inc);
                pji_u += bsz2;
                bsz_cols += bsz_col;
            }

            if (nrj_u > 1) {
                SGEMV(&chn, &bsz_row, &bsz_cols, &minusone, &sys.Mg_sp[pj_u], &bsz_row, &rDenseRow_sp[0],
                      &inc, &one, &r_sp[pjjR], &inc);
            }
            
            // Compute r_j <- U_jj \ r_j
            SGEMV(&chn, &bsz_row, &bsz_row, &one, &sys.Mg_sp[pji_u], &bsz_row, &r_sp[pjjR],
                  &inc, &zero, &sys.rtmp_sp[0], &inc);
            SCOPY(&bsz_row, &sys.rtmp_sp[0], &inc, &r_sp[pjjR], &inc);
        }
    }
    
//     for (Int i = 0; i < numIntRows*bsz; i++)
//         r[i] = (double) r_sp[i];
    sys.r.insert(sys.r.begin(), sys.r_sp.begin(), sys.r_sp[blkStartR[numIntRows+1]]);
}

void applyBILUk_dp(sysstruct &sys, double * r, Int * requestCounter)
{
    // Mg is stored in memory in the same order as is required for the BILUk solve (L first, then U).

    USE _sp ONCE _sp IS VALIDATED
    
#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && sys.parallelPreconditioner == ?? && (sys.reorderMethod == 0 || sys.reorderMethod == 1 || sys.reorderMethod == 2)) {
//     if (sys.nproc > 1 && preconditioner == 0 && (sys.reorderMethod == 0 || sys.reorderMethod == 1 || sys.reorderMethod == 2)) {
        recvvector(sys, &r[0], requestCounter);
    }
#endif
    
    if (sys.reorderMethod == 3 || sys.reorderMethod == 4) {
        error("Mixed-precision algorithm not implemented for MDF with constraints.\n");
    }
    
    Int inc = 1, info;
    Int nrowsP = sys.numRowsP;
    Int numIntRows = sys.numIntRows;
    double one = 1.0, zero = 0.0, minusone = -1.0;
    char chn = 'N';
    
    Int *ordered2unordered = &sys.ordered2unordered[0];
    Int *Uoent2entStart = &sys.Uoent2entStart[0];
    Int *Loent2entStart = &sys.Loent2entStart[0];
    Int *LUoent2ent = &sys.LUoent2ent[0];
    Int *UorowStart = &sys.UorowStart[0];
    Int *LorowStart = &sys.LorowStart[0];
    Int *rowStartR = &sys.rowStartR[0];
    Int *entSz = &sy.entSz[0];
    double *rDenseRow = &sys.xDenseRow[0][0];
    
    // Forward solve r <- L \ r
    if (sys.precSolveImplementation == 0) {
        for (Int rj=0; rj<numIntRows; rj++) {
            Int rjj = ordered2unordered[rj];
            Int pj_l = LorowStart[rj];
            Int pjSp_l = Loent2entStart[rj];
            Int nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];
            
            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            Int pji_l = pj_l;
            for (Int i=0; i<nrj_l; i++) {
                //pji_l = pj_l+i;
                Int rii = LUoent2ent[pjSp_l+i];
                Int piiR = rowStartR[rii];
                Int bsz_col = entSz[rii];
                Int bsz2 = bsz_row*bsz_col;
                
                DGEMV(&chn, &bsz_row, &bsz_col, &minusone, &sys.Mg[pji_l], &bsz_row, &r[piiR],
                      &inc, &one, &r[pjjR], &inc);
                pji_l += bsz2;
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (Int rj=0; rj<numIntRows; rj++) {
            Int rjj = ordered2unordered[rj];
            Int pj_l = LorowStart[rj];
            Int pjSp_l = Loent2entStart[rj];
            Int nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            Int pji_l = pj_l;
            Int bsz_cols = 0;
            for (Int i=0; i<nrj_l; i++) {
                Int pjiSp_l = pjSp_l+i;
                Int rii = LUoent2ent[pjiSp_l];
                Int bsz_col = entSz[rii];
                Int piiR = rowStartR[rii];
                
                for (Int jjj = 0; jjj < bsz_col; jjj++)
                   rDenseRow[bsz_cols+jjj] = r[piiR+jjj];
                //DCOPY(&bsz, &r[piiR], &inc, &rDenseRow[bsz_cols], &inc);
                pji_l += bsz2;
                bsz_cols += bsz_col;
            }
            
            if (nrj_l > 0) {
                DGEMV(&chn, &bsz_row, &bsz_cols, &minusone, &sys.Mg[pj_l], &bsz_rpw, &rDenseRow[0],
                      &inc, &one, &r[pjjR], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }

#ifdef  HAVE_MPI
    // Receive input vector from entities not in the subdomain
    if (sys.nproc > 1 && parallelPreconditioner == ?? && (sys.reorderMethod == 3 || sys.reorderMethod == 4)) {
//     if (sys.nproc > 1 && preconditioner == 0 && (sys.reorderMethod == 0 || sys.reorderMethod == 1 || sys.reorderMethod == 2)) {
        recvvector(sys, &r[0], requestCounter);
    }
#endif
    
    if (sys.precSolveImplementation == 0) {
        for (Int rj=numIntRows; rj<nrowsP; rj++) {
            Int rjj = ordered2unordered[rj];
            Int pj_l = LorowStart[rj];
            Int pjSp_l = Loent2entStart[rj];
            Int nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];
            
            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            Int pji_l = pj_l;
            for (Int i=0; i<nrj_l; i++) {
                //pji_l = pj_l+i;
                Int rii = LUoent2ent[pjSp_l+i];
                Int piiR = rowStartR[rii];
                Int bsz_col = entSz[rii];
                Int bsz2 = bsz_row*bsz_col;
                
                DGEMV(&chn, &bsz_row, &bsz_col, &minusone, &sys.Mg[pji_l], &bsz_row, &r[piiR],
                      &inc, &one, &r[pjjR], &inc);
                pji_l += bsz2;
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (Int rj=numIntRows; rj<nrowsP; rj++) {
            Int rjj = ordered2unordered[rj];
            Int pj_l = LorowStart[rj];
            Int pjSp_l = Loent2entStart[rj];
            Int nrj_l = Loent2entStart[rj+1]-Loent2entStart[rj];
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];

            // Compute r_j <- r_j - L_ji*r_i for all i<j neighboring j
            Int pji_l = pj_l;
            Int bsz_cols = 0;
            for (Int i=0; i<nrj_l; i++) {
                Int pjiSp_l = pjSp_l+i;
                Int rii = LUoent2ent[pjiSp_l];
                Int bsz_col = entSz[rii];
                Int piiR = rowStartR[rii];
                
                for (Int jjj = 0; jjj < bsz_col; jjj++)
                   rDenseRow[bsz_cols+jjj] = r[piiR+jjj];
                //DCOPY(&bsz, &r[piiR], &inc, &rDenseRow[bsz_cols], &inc);
                pji_l += bsz2;
                bsz_cols += bsz_col;
            }
            
            if (nrj_l > 0) {
                DGEMV(&chn, &bsz_row, &bsz_cols, &minusone, &sys.Mg[pj_l], &bsz_rpw, &rDenseRow[0],
                      &inc, &one, &r[pjjR], &inc);
            }
            // r_j <- L_jj \ r_j does not need to be performed since the diagonal blocks of L are identity blocks
        }
    }
                    
    // Backward solve r <- U \ r
    if (sys.precSolveImplementation == 0) {
        for (Int j=0; j<nrowsP; j++) {
            Int rj = nrowsP-j-1;
            Int rjj = ordered2unordered[rj];
            Int pjSp = Uoent2entStart[rj]-1;
            Int pjSp_u = Uoent2entStart[rj+1];
            Int pj_u = UorowStart[rj+1];
            Int nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];
            
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];
            
            Int pji_u = pj_u;
            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            for (Int i=0; i<(nrj_u-1); i++) {
                Int rii = LUoent2ent[pjSp_u+i];
                Int piiR = rowStartR[rii];
                Int bsz_col = entSz[rii];
                Int bsz2 = bsz_row * bsz_col;
                
                DGEMV(&chn, &bsz_row, &bsz_col, &minusone, &sys.Mg[pji_u], &bsz_row, &r[piiR],
                      &inc, &one, &r[pjjR], &inc);
                pji_u += bsz2;
            }
            
            // Compute r_j <- U_jj \ r_j
            DGEMV(&chn, &bsz_row, &bsz_row, &one, &sys.Mg[pji_u], &bsz_row, &r[pjjR],
                  &inc, &zero, &sys.rtmp[0], &inc);
            
            DCOPY(&bsz_row, &sys.rtmp[0], &inc, &r[pjjR], &inc);
        }
    }
    else if (sys.precSolveImplementation == 1) {
        for (Int j=0; j<nrowsP; j++) {
            Int rj = nrowsP-j-1;
            Int rjj = ordered2unordered[rj];
            Int pjSp = Uoent2entStart[rj]-1;
            Int pjSp_u = Uoent2entStart[rj+1];
            Int pj_u = UorowStart[rj+1];
            Int nrj_u = Uoent2entStart[rj]-Uoent2entStart[rj+1];
            
            Int pjjR = rowStartR[rjj];
            Int bsz_row = entSz[rjj];

            Int pji_u = pj_u;
            Int bsz_cols = 0;
            // Compute r_j <- r_j - U_ji*r_i for all i>j neighboring j
            for (Int i=0; i<(nrj_u-1); i++) {
                Int rii = LUoent2ent[pjSp_u+i];
                Int piiR = rowStartR[rii];
                Int bsz_col = entSz[rii];
                Int bsz2 = bsz_row * bsz_col;
                Int pjiSp_u = pjSp_u+i;
                
                for (Int jjj = 0; jjj < bsz_col; jjj++)
                    rDenseRow[bsz_cols+jjj] = r[piiR+jjj];
                //SCOPY(&bsz_col, &r[piiR], &inc, &rDenseRow[bsz_cols], &inc);
                pji_u += bsz2;
                bsz_cols += bsz_col;
            }

            if (nrj_u > 1) {
                Int nb = bsz*(nrj_u-1);
                DGEMV(&chn, &bsz_row, &bsz_cols, &minusone, &sys.Mg[pj_u], &bsz_row, &rDenseRow[0],
                      &inc, &one, &r[pjjR], &inc);
            }

            // Compute r_j <- U_jj \ r_j
            DGEMV(&chn, &bsz_row, &bsz_row, &one, &sys.Mg[pji_u], &bsz_row, &r[pjjR],
                  &inc, &zero, &sys.rtmp[0], &inc);
            DCOPY(&bsz, &sys.rtmp[0], &inc, &r[pjjR], &inc);
        }
    }
}

void applyEntityInvBJ_sp(sysstruct &sys, float *r_sp, double *r)
{
    Int inc = 1;
    float one = 1.0, zero = 0.0, minusone = -1.0;
    char chn = 'N';
    
    Int nrowsP = sys.numRowsP;
    Int *rowStartP = &sys.rowStartP[0];
    Int *rowStartR = &sys.rowStartR[0];
    
// // // //     #pragma omp parallel for num_threads(sys.noThreads)
// // // //     for (Int i = 0; i < nrowsP*bsz; i++)
// // // //         r_sp[i] = (float) r[i];
    sys.r_sp.insert(sys.r_sp.begin(), sys.r.begin(), sys.r.end());
    
//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for private(pj)
        for (Int rj=0; rj<nrows; rj++) {
            Int pjP = rowStartJ[rj];
            Int pjR = rowStartR[rj];
            
            Int bsz = sys.entSz[rj];
            Int bsz2 = bsz*bsz;
            
            SGEMV(&chn, &bsz, &bsz, &one, &sys.Mg_sp[pjP], &bsz, &r_sp[pjR],
                  &inc, &zero, &sys.rtmps_sp[this_thread][0], &inc);
            SCOPY(&bsz, &sys.rtmps_sp[this_thread][0], &inc, &r_sp[pjR], &inc);
        }
//     }
    
// // // //     #pragma omp parallel for num_threads(sys.noThreads)
// // // //     for (Int i = 0; i < nrowsP*bsz; i++)
// // // //         r[i] = (double) r_sp[i];
    sys.r.insert(sys.r.begin(), sys.r_sp.begin(), sys.r_sp.end());
}

void applyEntityInvBJ_dp(sysstruct &sys, double *r)
{
    Int inc = 1;
    double one = 1.0, zero = 0.0, minusone = -1.0;
    char chn = 'N';
    
    Int nrowsP = sys.numRowsP;
    Int *rowStartP = &sys.rowStartP[0];
    Int *rowStartR = &sys.rowStartR[0];

//     #pragma omp parallel num_threads(sys.noThreads)
//     {
//         int this_thread = omp_get_thread_num();
        int this_thread = 0;
//         #pragma omp for private(pj)
        for (Int rj=0; rj<nrowsP; rj++) {
            Int pjP = rowStartJ[rj];
            Int pjR = rowStartR[rj];
            
            Int bsz = sys.entSz[rj];
            Int bsz2 = bsz*bsz;
    
            DGEMV(&chn, &bsz, &bsz, &one, &sys.Mg[pjP], &bsz, &r[pjR],
                  &inc, &zero, &sys.rtmps[this_thread][0], &inc);
            DCOPY(&bsz, &sys.rtmps[this_thread][0], &inc, &r[pjR], &inc);
        }
//     }
}

void computePreconditioner(sysstruct &sys, Int subdomainPreconditioner)
{
    if (subdomainPreconditioner == -1) {                 // No preconditioner
    }
    else if (subdomainPreconditioner == 0)      // BILU(k)
        computeBILUk(sys);
    else if (subdomainPreconditioner == 1)      // Entity-based block Jacobi
        computeEntityInvBJ(sys);
    
    // Convert preconditioner into single precision, if necessary:
    if (subdomainPreconditioner != -1 && sys.precPrecision == 0) {
//         for (Int i = 0; i < sys.Mg.size(); i++)
//             sys.Mg_sp[i] = (float) sys.Mg[i];
        sys.Mg_sp.insert(sys.Mg_sp.begin(), sys.Mg.begin(), sys.Mg.end());
    }
}

void applyPreconditioner_sp(sysstruct &sys, float *r_sp, double *r, Int subdomainPreconditioner, Int *requestCounter)
{   
    if (subdomainPreconditioner == -1) {        // No preconditioner
    }
    else if (subdomainPreconditioner == 0)      // BILU(k)
        applyBILUk_sp(sys, r_sp, r, requestCounter);
    else if (subdomainPreconditioner == 1)      // Entity-based block Jacobi
        applyEntityInvBJ_sp(sys, r_sp, r);
}

void applyPreconditioner_dp(sysstruct &sys, double *r, Int subdomainPreconditioner, Int *requestCounter)
{   
    if (subdomainPreconditioner == -1) {        // No preconditioner
    }
    else if (subdomainPreconditioner == 0)      // BILU(k)
        applyBILUk_dp(sys, r, requestCounter);
    else if (subdomainPreconditioner == 1)      // Entity-based block Jacobi
        applyEntityInvBJ_dp(sys, r);
}

void applyPreconditioner(sysstruct &sys, double *r, Int subdomainPreconditioner, Int *requestCounter)
{   
    if (sys.precPrecision == 0 && sys.robustMode == 0) {
        float *r_sp = &sys.r_sp[0];
        applyPreconditioner_sp(sys, r_sp, r, subdomainPreconditioner, requestCounter);
    }
    else
        applyPreconditioner_dp(sys, r, subdomainPreconditioner, requestCounter);
}

void computeEnt2entP_BILUk(sys)
{
    // Compute ent2entP and ent2entStartP for BILU(k) preconditioner
    computeKlevelConnectivity(&sys.ent2entP, &sys.ent2entStartP, &sys.ent2entJ, &sys.ent2entStartJ, sys.BILUkLevel);
}

void computeEnt2entP_EntityInvBJ(sys)
{
    // Compute ent2entP and ent2entStartP for entity-based block Jacobi preconditioner
    
    Int numIntEnt = sys.numIntEnt;
    
    sys.ent2entP.resize(numIntEnt);
    sys.ent2entStartP.resize(numIntEnt+1);
    for (Int i=0; i<numIntEnt; i++) {
        sys.ent2entP[i] = i;
        sys.ent2entStartP[i] = i;
    }
    sys.ent2entStartP[numIntEnt] = numIntEnt;
}

void computeEnt2entP_NoPrec(sys)
{
    // Compute ent2entP and ent2entStartP for no preconditioner. The sparsity pattern of the identity matrix is used
    
    Int numIntEnt = sys.numIntEnt;
    
    sys.ent2entP.resize(numIntEnt);
    sys.ent2entStartP.resize(numIntEnt+1);
    for (Int i=0; i<numIntEnt; i++) {
        sys.ent2entP[i] = i;
        sys.ent2entStartP[i] = i;
    }
    sys.ent2entStartP[numIntEnt] = numIntEnt;
}

void computeEnt2entP(sys)
{
    // Compute ent2entP and ent2entStartP for subdomain preconditioner
    
    if (sys.subdomainPreconditioner == -1)          // No preconditioner
        computeEnt2entP_NoPrec(sys);
    else if (sys.subdomainPreconditioner == 0)      // BILU(k)
        computeEnt2entP_BILUk(sys);
    else if (sys.subdomainPreconditioner == 1)      // Entity-based block Jacobi
        computeEnt2entP_EntityInvBJ(sys)
}

void preprocessSparsityPreconditioner(sys)
{
    computeEnt2entP(sys);
    
    if (sys.subdomainPreconditioner == 0)
        computeMgMinusHg(sys);
}

#endif
