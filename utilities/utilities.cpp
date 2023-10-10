#ifndef __UTILITIES
#define __UTILITIES

// void splitArray1D(Int arrayLen, Int this_thread, Int noThreads, Int* chunkLen, Int* startPosition)
// {
//     *startPosition = ((arrayLen * this_thread) / noThreads);
//     Int endPosition = ((arrayLen * (this_thread+1)) / noThreads);
//     
//     if (this_thread == (noThreads-1)) {   // This is reduntant unless we go over the range in which double represent integers exactly
//         if (endPosition != arrayLen) {
//             printf("Warning 4TG67N in splitArray1D.\n");
//             endPosition = arrayLen;
//         }
//     }
//     
//     *chunkLen = endPosition - *startPosition;
// }

Int getElemtype(Int nd, Int nve)
{
    Int elemtype;
    
    if (nd  == 1 && nve == 2)       // Line
        elemtype = 0;
    else if (nd == 2 && nve == 3) // Tri
        elemtype = 0;
    else if (nd == 2 && nve == 4) // Quad
        elemtype = 1;
    else if (nd == 3 && nve == 4)	// Tet
        elemtype = 0;
    else if (nd == 3 && nve == 8)	// Hex
        elemtype = 1;
    else
        error("The number of vertices per element does not match any implemented element type.\n");

    return elemtype;
}

Int getSpaceDimension(Int porder, Int nd, Int elemtype)
{
    Int dim, i;
    if (elemtype == 0) {
        Int num = 1, den = 1;
        for (i = 0; i < nd; i++)
            num *= (porder + i + 1);
        for (i = 0; i < nd; i++)
            den *= (i + 1);
        if ((num % den) != 0)
            error("Error RF5T78 in getSpaceDimension.\n");
        else
            dim = num / den;
    }
    else if (elemtype == 1) {
        dim = 1;
        for (Int i = 0; i < nd; i++)
            dim *= porder+1;
    }
    else
        error("elemtype has invalid value.");
    
    return dim;
}

void getNodalIndexOfVertex(Int *vertexIndex, Int vertexNo, Int nd, Int elemtype, Int porder)
{
    if (vertexNo == 0)
        *vertexIndex = 1 - 1;
    else if (vertexNo == 1)
        *vertexIndex = porder + 1 - 1;
    else if (vertexNo == 2) {
        if (elemtype == 0)
            *vertexIndex = ((porder+1)*(porder+2))/2 - 1;
        else if (elemtype == 1)
            *vertexIndex = (porder+1)*(porder+1) - 1;
    }
    else if (vertexNo == 3) {
        if (elemtype == 0)
            *vertexIndex = ((porder+1)*(porder+2)*(porder+3))/6 - 1;
        else if (elemtype == 1)
            *vertexIndex = (porder+1)*porder + 1 - 1;
    }
    else if (vertexNo == 4) {
        if (elemtype == 0)
            error("Error GB7UJ6 in getNodalIndexOfVertex.\n");
        else if (elemtype == 1)
            *vertexIndex = (porder+1)*(porder+1)*porder + 1 - 1;
    }
    else if (vertexNo == 5) {
        if (elemtype == 0)
            error("Error GB7UJ7 in getNodalIndexOfVertex.\n");
        else if (elemtype == 1)
            *vertexIndex = (porder+1)*(porder+1)*porder + porder + 1 - 1;
    }
    else if (vertexNo == 6) {
        if (elemtype == 0)
            error("Error GB7UJ8 in getNodalIndexOfVertex.\n");
        else if (elemtype == 1)
            *vertexIndex = (porder+1)*(porder+1)*(porder+1) - 1;
    }
    else if (vertexNo == 7) {
        if (elemtype == 0)
            error("Error GB7UJ9 in getNodalIndexOfVertex.\n");
        else if (elemtype == 1)
            *vertexIndex = (porder+1)*(porder+1)*porder + (porder+1)*porder + 1 - 1;
    }
}

void getProjectionMatrix(Int *ndims, Int projporder, double *P)
{
    // Compute projection matrix to a low-p subspace
    // Note: The orthogonal basis are assumed to be ordered in non-decreasing polynomial order. (The DGEMM line below is coded based on that assumption.)
    
    Int nd = ndims[0];
    Int nve = ndims[3];
    Int npv = ndims[9];
    Int porder = ndims[15];
    
    Int N = npv * npv;
    Int elemtype = getElemtype(nd, nve); 
    Int projDim = getSpaceDimension(projporder, nd, elemtype);      // Get dimensionality of projection subspace
    
    double zero = 0.0, one = 1.0;
    char chn = 'N';
    string filename;
    
    if (projporder >= porder) {     // P is just the identity matrix
        for (Int i = 0; i < npv*npv; i++)
            P[i] = 0.0;
        for (Int i = 0; i < npv; i++)
            P[i*(npv+1)] = 1.0;
    }
    else {
        vector<double> orth2nodal;
        vector<double> nodal2orth;
        orth2nodal.resize(N);
        nodal2orth.resize(N);
        
        // Get filename
        if (elemtype == 0) {
            if (nd == 1)
                filename = "basisChange_line_p" + NumberToString(porder) + ".bin";
            else if (nd == 2)
                filename = "basisChange_tri_p" + NumberToString(porder) + ".bin";
            else if (nd == 3)
                filename = "basisChange_tet_p" + NumberToString(porder) + ".bin";
        }
        else if (elemtype == 1) {
            if (nd == 1)
                filename = "basisChange_line_p" + NumberToString(porder) + ".bin";
            else if (nd == 2)
                filename = "basisChange_quad_p" + NumberToString(porder) + ".bin";
            else if (nd == 3)
                filename = "basisChange_hex_p" + NumberToString(porder) + ".bin";
        }
        else
            error("elemtype not recognized.\n");
        
        // Read change of basis matrices
        ifstream in(filename.c_str(), ios::in | ios::binary);
        if (!in) {
            error("Unable to open file " + filename);
        }
        else {
            readdarray(in, orth2nodal, N);
            readdarray(in, nodal2orth, N);
        }
        in.close();
        
        // Compute projection matrix:
        DGEMM(&chn, &chn, &npv, &npv, &projDim, &one, &orth2nodal[0], &npv, &nodal2orth[0], &npv, &zero, &P[0], &npv);
    }
}

void projectUDGfields(double *projectedFields, double *P, double *UDG, Int *fields2Project, Int numFields2Proj, Int npv)
{
    Int inc = 1, initUDG;
    double one = 1.0, zero = 0.0;
    char chn = 'N';
    
    for (Int i = 0; i < numFields2Proj; i++) {
        initUDG = fields2Project[i] * npv;
        DGEMV(&chn, &npv, &npv, &one, &P[0], &npv, &UDG[initUDG], &inc, &zero, &projectedFields[i*npv], &inc);
    }
}

void projectFields(double *projectedFields, double *fields, double *projMatrix, Int numFields2Proj, Int npv)
{
    Int inc = 1;
    double one = 1.0, zero = 0.0;
    char chn = 'N';
    
    for (Int i = 0; i < numFields2Proj; i++)
        DGEMV(&chn, &npv, &npv, &one, &projMatrix[0], &npv, &fields[i*npv], &inc, &zero, &projectedFields[i*npv], &inc);
}

void projectField(double *projectedField, double *field, double *projMatrix, Int npv)
{
    Int inc = 1;
    double one = 1.0, zero = 0.0;
    char chn = 'N';
    
    DGEMV(&chn, &npv, &npv, &one, &projMatrix[0], &npv, &field[0], &inc, &zero, &projectedField[0], &inc);
}

#endif
