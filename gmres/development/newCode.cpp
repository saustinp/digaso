
// Written by: C. Nguyen & P. Fernandez

void block_crs_k(vector <Int> *colind_K_p , vector <Int> *rowpts_K_p , 
                 vector <Int> *colind_p   , vector <Int> *rowpts_p   , Int Klevel)
{
    // Compute block_crs format for K-level overlap.
    // This is necessary to compute the sparsity pattern of BILU(k) preconditioners.
    
    // Note:
    // 1) entity (in hybridized DG context) = node (in graph theory context)
    // 2) The input and output pointers CAN be pointing to the same vector.
    
    if (Klevel < 0)
        error("block_crs_k only valid for non-negative level overlaps.\n");
    
    std::vector<Int>::iterator it;
    
    Int nrows = rowpts_p[0].size() - 1;
    Int nedges_k = colind_p[0].size();
    
    vector <Int> rowpts_kMinus1(nrows+1,0);
    vector <Int> colind_kMinus1;
    vector <Int> colind_i_k;
    
    // (colind, rowpts) -> (colind_k, rowpts_k):
    vector <Int> rowpts_k(nrows+1, 0);
    for (int i=0; i<nrows+1; i++)
        rowpts_k[i] = rowpts_p[0][i];
    vector <Int> colind_k(nedges_k, 0);
    for (int i=0; i<nedges_k; i++)
        colind_k[i] = colind_p[0][i];
    
    Int maxEdgesPerEnt_k = 0;
    for (int i=0; i<nrows; i++)
        maxEdgesPerEnt_k = max(maxEdgesPerEnt_k, rowpts_k[i+1] - rowpts_k[i]);
    
    // Recursively compute K-level connectivity graph:
    for (int k = 0; k < Klevel; k++) {
        // (colind_k, rowpts_k) -> (colind_kMinus1, rowpts_kMinus1):
        Int nedges_kMinus1 = nedges_k;
        colind_kMinus1.resize(nedges_kMinus1);
        for (int i=0; i<nedges_kMinus1; i++)
            colind_kMinus1[i] = colind_k[i];
        for (int i=0; i<nrows+1; i++)
            rowpts_kMinus1[i] = rowpts_k[i];
        Int maxEdgesPerEnt_kMinus1 = maxEdgesPerEnt_k;
        
        // Allocate enough memory for k-level connectivity vector:
        colind_k.resize( ((maxEdgesPerEnt_kMinus1 / 4) + 1) * nedges_kMinus1);
        
        // Loop over all entities:
        nedges_k = 0;
        rowpts_k[0] = 0;
        colind_i_k.resize(maxEdgesPerEnt_kMinus1 * maxEdgesPerEnt_kMinus1);
        maxEdgesPerEnt_k = 0;
        for (int i = 0; i < nrows; i++) {
            Int nedges_i_kMinus1 = rowpts_kMinus1[i+1] - rowpts_kMinus1[i];
            Int start_i = rowpts_kMinus1[i];
            
            // Compute colind for i-th entity in k-level connectivity graph:
            it = colind_i_k.begin();
            for (int jj = 0; jj < nedges_i_kMinus1; jj++) {
                Int j = colind_kMinus1[start_i+jj];
                Int start_j = rowpts_kMinus1[j];
                Int end_j = rowpts_kMinus1[j+1];
                it = colind_i_k.insert(it, &colind_kMinus1[start_j], &colind_kMinus1[end_j]);
            }
            
            std::sort (colind_i_k.begin(), it);
            it = std::unique (colind_i_k.begin(), it);
            Int nedges_i_k = it-colind_i_k.begin();
            if (nedges_i_k < 1)
                error("nedges_i_k < 1 in block_crs_k was detected.\n");
            maxEdgesPerEnt_k = max(maxEdgesPerEnt_k, nedges_i_k);
            
            // Update colind_k and rowpts_k with connectivities of i-th entity:
            nedges_k += nedges_i_k;
            rowpts_k[i+1] = rowpts_k[i] + nedges_i_k;
            if (colind_k.size() < nedges_k)
                colind_k.resize(nedges_k);
            for (int j=0; j < nedges_i_k; j++)
                colind_k[rowpts_k[i] + j] = colind_i_k[j];
        }
        colind_k.resize(nedges_k);
    }
    
    // (colind_k, rowpts_k) -> (colind_K, rowpts_K):
    colind_K_p[0].resize(nedges_k);
    for (int i=0; i<nedges_k; i++)
        colind_K_p[0][i] = colind_k[i];
    rowpts_K_p[0].resize(nrows+1);
    for (int i=0; i<nrows+1; i++)
        rowpts_K_p[0][i] = rowpts_k[i];
    
//     // If colind_K_p is a pointer to a vector of vectors (new desired implementation)...
//     for (int i=0; i<nrows; i++)
//         copy colind_k[rowpts_k[i]] to colind_K_p[0][i]
}



Int matchSparsity(sysstruct &sys, Int ri, Int j)
{
    Int piP = sys.rowptsP[ri];
    Int piJ = sys.ent2entStart[ri];
    Int *rowptsJ = &sys.ent2entStart[0];
    Int rj = sys.colindP[piP+j];
    Int *colindJ = &sys.ent2ent[0];
    
    Int match = 0;
    Int nriJ = rowptsJ[ri+1] - rowptsJ[ri];
    for (Int k = 0; k < nriJ; k++) {
        if (colindJ[piJ+k] == rj) {
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
    
// // //     if (sys.subdomainPreconditioner != 0)
// // //         error("computeMgMinusHg only valid for BILU(k) subdomain preconditioner.\n");
    
    Int *rowptsP = &sys.rowptsP[0];
    
    sys.MgMinusHg.resize(sys.colindP.size());
    
    Int nrowsP = sys.nrowsP;
    Int count = 0;
    for (Int ri = 0; ri < nrowsP; ri++) {
        Int nriP = rowptsP[ri+1] - rowptsP[ri];
        for (Int i = 0; i < nriP; i++)
            sys.MgMinusHg[count++] = matchSparsity(sys, ri, i);
    }
    
    if (count != sys.colindP.size())
        error("Error H4VCB0 in computeMgMinusHg.\n");
}

void copyHg2MgWithMgFormat(sysstruct &sys)
{
    Int inc = 1;
    Int nrowsP = sys.nrowsP;
    
    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    
    Int *rowptsJ = &sys.ent2entStart[0];
    Int *rowptsP = &sys.rowptsP[0];
    Int *colindP = &sys.colindP[0];
    Int *Lorowpts = &sys.Lorowpts[0];
    Int *Uorowpts = &sys.Uorowpts[0];
    
    for (Int rjj=0; rjj<nrowsP; rjj++) {
        Int rj = sys.unordered2ordered[rjj];
        
        Int pjP   = sys.Uoent2entStart[rj]-1;
        Int pjP_l = sys.Loent2entStart[rj];
        Int pjP_u = sys.Uoent2entStart[rj+1];
        
        Int pjjJ  = rowptsJ[rjj];
        Int pjjP  = rowptsP[rjj];
        
        Int nrjJ  = rowptsJ[rjj+1]-rowptsJ[rjj];
        Int nrjP  = rowptsP[rjj+1]-rowptsP[rjj];
        Int nrj_l = sys.Loent2entStart[rj+1]-sys.Loent2entStart[rj];        // Neighbors in L
        Int nrj_u = sys.Uoent2entStart[rj]-sys.Uoent2entStart[rj+1];        // Neighbors in U (includes the diagonal block)
        
        // Copy diagonal block:
        DCOPY(&bsz2, &sys.Hg[bsz2*pjjJ], &inc, &sys.Mg[bsz2*pjP], &inc);
        
        // Copy off-diagonal blocks:
        Int i_u = 0, i_l = 0;
        Int pjiP_l = pjP_l;
        Int pjiP_u = pjP_u;
        Int pjjiJ = pjjJ + 1;
        for (Int i=1; i<nrjP; i++) {
            rii = colindP[pjjP+i];
            ri = sys.unordered2ordered[rii];
            
            if (ri > rj) {
                // U block
                i_u ++;
                if (! sys.MgMinusHg[rowptsP[rjj]+i])
                    std::fill(&sys.Mg[bsz2*pjiP_u], &sys.Mg[bsz2*(pjiP_u+1)], 0.0);
                else {
                    DCOPY(&bsz2, &sys.Hg[bsz2*pjjiJ], &inc, &sys.Mg[bsz2*pjiP_u], &inc);
                    pjjiJ += 1;
                }
                pji_u += 1;
            }
            else if (ri < rj) {
                // L block
                i_l ++;
                if (! sys.MgMinusHg[rowptsP[rjj]+i])
                    std::fill(&sys.Mg[bsz2*pjiP_l], &sys.Mg[bsz2*(pjiP_l+1)], 0.0);
                else {
                    DCOPY(&bsz2, &sys.Hg[bsz2*pjjiJ], &inc, &sys.Mg[bsz2*pjiP_l], &inc);
                    pjjiJ += 1;
                }
                pji_l += 1;
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

