#ifndef __BLOCK_CRS_K
#define __BLOCK_CRS_K

// Written by: C. Nguyen & P. Fernandez

void block_crs_k(vector <Int> *colind_K_p , vector <Int> *rowpts_K_p , 
                 vector <Int> *colind_p   , vector <Int> *rowpts_p   , Int maxRow, Int Klevel)
{
    // Compute block_crs format for K-level overlap.
    // This is necessary to compute the sparsity pattern of BILU(k) preconditioners.
    
    // Note:
    // 1) entity (in hybridized DG context) = node (in graph theory context)
    // 2) The input and output pointers CAN be pointing to the same vector.
    
    if (Klevel < 0)
        error("block_crs_k only valid for non-negative level overlaps.\n");
    if (maxRow > rowpts_p[0].size() - 1)
        error("maxRow > nrows detected in block_crs_k.\n");
    
    Int nedges_k = colind_p[0].size();
    
    vector <Int> rowpts_kMinus1(maxRow+1,0);
    vector <Int> colind_kMinus1;
    vector <Int> colind_i_k;
    vector <Int> colind_i_k_tmp;
    
    // (colind, rowpts) -> (colind_k, rowpts_k):
    vector <Int> rowpts_k(maxRow+1, 0);
    vector <Int> colind_k(nedges_k, 0);
    nedges_k = 0;
    rowpts_k[0] = 0;
    for (int i=0; i<maxRow; i++) {
        Int nri = rowpts_p[0][i+1] - rowpts_p[0][i];
        Int nedges_i_k = 0;
        for (int j=0; j<nri; j++) {
            Int rj = colind_p[0][rowpts_p[0][i] + j];
            if(rj < maxRow) {
                colind_k[nedges_k] = rj;
                nedges_k++;
                nedges_i_k++;
            }
        }
        rowpts_k[i+1] = rowpts_k[i] + nedges_i_k;
    }
    colind_k.resize(nedges_k);
    
    Int maxEdgesPerEnt_k = 0;
    for (int i=0; i<maxRow; i++)
        maxEdgesPerEnt_k = max(maxEdgesPerEnt_k, rowpts_k[i+1] - rowpts_k[i]);
    
    // Recursively compute K-level connectivity graph:
    for (int k = 0; k < Klevel; k++) {
        // (colind_k, rowpts_k) -> (colind_kMinus1, rowpts_kMinus1):
        Int nedges_kMinus1 = nedges_k;
        colind_kMinus1.resize(nedges_kMinus1);
        for (int i=0; i<nedges_kMinus1; i++)
            colind_kMinus1[i] = colind_k[i];
        for (int i=0; i<maxRow+1; i++)
            rowpts_kMinus1[i] = rowpts_k[i];
        Int maxEdgesPerEnt_kMinus1 = maxEdgesPerEnt_k;
        
        // Allocate enough memory for k-level connectivity vector:
        colind_k.resize( 10 * nedges_kMinus1);
//         colind_k.resize( ((maxEdgesPerEnt_kMinus1 / 4) + 1) * nedges_kMinus1);
        
        // Loop over all entities:
        nedges_k = 0;
        rowpts_k[0] = 0;
        colind_i_k.resize(maxEdgesPerEnt_kMinus1 * maxEdgesPerEnt_kMinus1);
        maxEdgesPerEnt_k = 0;
        for (int i = 0; i < maxRow; i++) {
            Int nedges_i_kMinus1 = rowpts_kMinus1[i+1] - rowpts_kMinus1[i];
            Int start_i = rowpts_kMinus1[i];
            
            // Compute colind for i-th entity in k-level connectivity graph:
            Int count = 0;
            for (int jj = 0; jj < nedges_i_kMinus1; jj++) {
                Int j = colind_kMinus1[start_i+jj];
                Int start_j = rowpts_kMinus1[j];
                Int end_j = rowpts_kMinus1[j+1];
                colind_i_k.insert(colind_i_k.begin()+count, &colind_kMinus1[start_j], &colind_kMinus1[end_j]);  // HERE
                count += (end_j-start_j);
// // //                 set_union(&colind_kMinus1[start_j],&colind_kMinus1[end_j],&colind_i_k[0],&colind_i_k[count],back_inserter(v3));
            }
            if (count < 1)
                error("Error No. 1 in block_crs_k.\n");
            
// // // //             Int nedges_i_k;
// // // //             uniqueiarray(&nedges_i_k, colind_i_k_tmp, colind_i_k, count);
// // // //             printf("%d, %d.\n", nedges_i_k, colind_i_k_tmp.size());
            uniqueiarray(colind_i_k_tmp,colind_i_k,count);      // HERE
            Int nedges_i_k = colind_i_k_tmp.size();
            if (nedges_i_k < 1)
                error("nedges_i_k < 1 in block_crs_k was detected.\n");
            
            colind_i_k[0] = i;
            Int indx = 1;
            for (int j=0; j<nedges_i_k; j++) {
                if (colind_i_k_tmp[j] != i) {
                    colind_i_k[indx] = colind_i_k_tmp[j];
                    indx++;
                }
            }
            if (indx != nedges_i_k)
                error("Error No. 1 in block_crs_k.\n");
            maxEdgesPerEnt_k = max(maxEdgesPerEnt_k, nedges_i_k);
            
            // Update colind_k and rowpts_k with connectivities of i-th entity:
            nedges_k += nedges_i_k;
            rowpts_k[i+1] = rowpts_k[i] + nedges_i_k;
            if (colind_k.size() < nedges_k)
                colind_k.resize( 2 * nedges_k );
            for (int j=0; j < nedges_i_k; j++)
                colind_k[rowpts_k[i] + j] = colind_i_k[j];
        }
        colind_k.resize(nedges_k);
    }
    
    // (colind_k, rowpts_k) -> (colind_K, rowpts_K):
    colind_K_p[0].resize(nedges_k);
    for (int i=0; i<nedges_k; i++)
        colind_K_p[0][i] = colind_k[i];
    rowpts_K_p[0].resize(maxRow+1);
    for (int i=0; i<maxRow+1; i++)
        rowpts_K_p[0][i] = rowpts_k[i];
}

#endif
