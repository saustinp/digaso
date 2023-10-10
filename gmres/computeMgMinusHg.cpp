#ifndef __COMPUTEMGMINUSHG
#define __COMPUTEMGMINUSHG

// Written by: C. Nguyen & P. Fernandez

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

void computeMgMinusHg(sysstruct &sys)
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

#endif
