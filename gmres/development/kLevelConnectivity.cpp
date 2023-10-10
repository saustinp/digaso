#ifndef __KLEVELCONNECTIVITY
#define __KLEVELCONNECTIVITY

// Written by: C. Nguyen & P. Fernandez

void block_crs_k(vector <Int> *colind_K_p, vector <Int> *rowpts_K_p, 
            vector <Int> *colind_p  , vector <Int> *rowpts_p  , Int Klevel)
{
    // Compute K-level connectivity graph from 0-level connecivity graph.
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

#endif
