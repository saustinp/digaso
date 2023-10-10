#ifndef __GETMESHT
#define __GETMESHT

// Written by: C. Nguyen & P. Fernandez

void get_mesh_t(meshstruct &mesh, appstruct &app, Int* ndims)
{
    // TODO: mesh.t is not correct if there are periodic boundaries
    
    Int i, j, k, ii, jj, ij, vertexIndex, found;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int nfe = ndims[2];
    Int nve = ndims[3];
    Int ne = ndims[5];
    Int nv = ndims[7];
    Int npv = ndims[9];
    Int npf = ndims[10];
    Int ndf = npf * nfe;
    Int porder = app.porder;
    Int elemtype = getElemtype(nd, nve);
    
    double tol = 1.0e-9;
    double dist = 0.0;
    
    double *pWithDuplications = new double [ne*nve*nd];
    mesh.t.resize(nve*ne);
    
    Int numVertices = 0;
    for (i = 0; i < ne; i++){
        for (j = 0; j < nve; j++) {
            getNodalIndexOfVertex(&vertexIndex, j, nd, elemtype, porder);
            for (k = 0; k < nd; k++)
                pWithDuplications[i*nve*nd+j*nd+k] = mesh.dgnodes[i*npv*ncd+k*npv+vertexIndex];
        }
    }
    for (i = 0; i < ne; i++) {
        for (j = 0; j < nve; j++) {
            found = 0;
            for (ij = 0; ij < i*nve; ij++) {
                ii = ij / nve;
                jj = ij % nve;
                dist = 0.0;
                for (k = 0; k < nd; k++)
                    dist += (pWithDuplications[i*nve*nd+j*nd+k]-pWithDuplications[ij*nd+k]) * (pWithDuplications[i*nve*nd+j*nd+k]-pWithDuplications[ij*nd+k]);
                dist = sqrt(dist);
                
                if (dist <= tol) {
                    mesh.t[i*nve+j] = mesh.t[ii*nve+jj];
                    found = 1;
                    break;
                }
            }
            if (found == 0) {
                mesh.t[i*nve+j] = numVertices;
                numVertices += 1;
            }
        }
    }

    delete[] pWithDuplications;

//    if (numVertices != nv) {
//        error("Number of vertices found not consistent with binary file input.\n");
//    }
}

#endif
