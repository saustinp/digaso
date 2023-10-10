#ifndef __DG2P1CG
#define __DG2P1CG

// Written by C. Nguyen and P. Fernandez

void DG_2_p1CG(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims, Int avgMethod)
{
    // avgMethod:   0: Vertex value given by David's reconstruction. 1: Vertex value is average of nodes on vertex at different elements
    
    Int nd = mesh.nd;
    Int nv = mesh.nv;
    Int ne = mesh.ne;
    Int porder = app.porder[??];
    double neighMearures;
    
    Int zeroPassiveElemets = 0;     // 1: Zero out elements for which the field is nonpositive on a nonempty portion of the element
    
    // Get number of total vertices in the mesh (vertices on periodic boundaries are counted twice, i.e., separately)
    Int numTotalVertices = 0;
    for (Int i = 0; i < mesh.t.size(); i++) {
        if (mesh.t[i] > numTotalVertices)
            numTotalVertices = mesh.t[i];
    }
    numTotalVertices += 1;
   if (nv != numTotalVertices)
       error("There are ghost vertices in the mesh.\n");

    // Creation of a matrix that contains the vertex-to-element connectivity:
    Int *numElementsAdj2Vertex = new Int[numTotalVertices];
    Int *vertex2elem = new Int[numTotalVertices*50];      // We assume at most 50 elements are connected to a vertex.
    Int *vertex2node = new Int[numTotalVertices*50];      // We assume at most 50 elements are connected to a vertex.
    double *fieldVertices = new double[numTotalVertices];
    Int *activeElement = new Int[ne];
    for (Int iv = 0; iv < numTotalVertices; iv++)
        numElementsAdj2Vertex[iv] = 0;
    for (Int ie = 0; ie < ne; ie++) {
        Int elemtype = mesh.elementype[ie];
        Int nve = mesh.elements[elemtype].nve;
        Int npv = mesh.elements[elemtype].npv;
                
        for (Int iv = 0; iv < nve; iv++) {
            Int currentVertex = mesh.t[iv*ne+ie];
            getNodalIndexOfVertex(&currentNode, iv, nd, elemtype, porder);
            numElementsAdj2Vertex[currentVertex] += 1;
            if (numElementsAdj2Vertex[currentVertex] > 50)
                error("Over 50 elements are connected to a vertex.\n");
            vertex2elem[currentVertex + numTotalVertices*(numElementsAdj2Vertex[currentVertex]-1)] = ie;
            vertex2node[currentVertex + numTotalVertices*(numElementsAdj2Vertex[currentVertex]-1)] = currentNode;
        }
        
        activeElement[ie] = 1;
        for (Int j = 0; j < npv; j++) {
            if (DGfield[sol.scalarFieldStart[ie] + j] <= 0.0) {
                activeElement[ie] = 0;
                break;
            }
        }
    }
    
    double *avgFieldElem = new double[ne];
    if (avgMethod == 0) {      // 0: Vertex value given by David's reconstruction
        // Compute average field in each element:
        for (Int ie = 0; ie < ne; ie++) {
            Int elemtype = mesh.elementype[ie];
            Int npv = mesh.elements[elemtype].npv;
            avgFieldElem[ie] = 0.0;
            for (Int j = 0; j < npv; j++) {
//                 for (Int k = 0; k < nqvR; k++)
//                     avgFieldElem[currentElem] += master.shapvgR[k*npv+j] * DGfield[currentElem*npv + j];
                avgFieldElem[ie] += DGfield[sol.scalarFieldStart[ie] + j];     // TODO: This is only an approximation to the true average. Some weighting is required for the exact average.
            }
            avgFieldElem[ie] /= ((double) npv);
        }
        
        // Compute field value at each vertex (average of the average field in the adjacent elements):
        for (Int iv = 0; iv < numTotalVertices; iv++) {
            fieldVertices[iv] = 0.0;
//             neighMearures = 0.0;
            for (Int j = 0; j < numElementsAdj2Vertex[iv]; j++) {
                Int currentElem = vertex2elem[iv + j*numTotalVertices];
                fieldVertices[iv] += avgFieldElem[currentElem];
//                 fieldVertices[currentVertex] += elemMeasure[currentElem] * avgFieldElem[currentElem];
//                 neighMearures += elemMeasure[currentElem];
            }
            fieldVertices[iv] /= max((double) numElementsAdj2Vertex[iv],1.0);
//             fieldVertices[currentVertex] /= max(neighMearures,1.0e-8);
        }
    }
    else if (avgMethod == 1) {     // 1: Vertex value is average of nodes on vertex at different elements
        // Compute field value at each vertex (average of the average field in the adjacent elements):
        for (Int iv = 0; iv < numTotalVertices; iv++) {
            fieldVertices[iv] = 0.0;
//             neighMearures = 0.0;
            for (Int j = 0; j < numElementsAdj2Vertex[iv]; j++) {
                Int currentElem = vertex2elem[iv + j*numTotalVertices];
                Int currentNode = vertex2node[iv + j*numTotalVertices];
                fieldVertices[iv] += DGfield[sol.scalarFieldStart[currentElem] + currentNode];
//                 fieldVertices[currentVertex] += elemMeasure[currentElem] * DGfield[currentElem*npv + currentNode];
//                 neighMearures += elemMeasure[currentElem];
            }
            fieldVertices[iv] /= max((double) numElementsAdj2Vertex[iv],1.0);
//             fieldVertices[currentVertex] /= max(neighMearures,1.0e-8);
        }
    }
    else
        error("Averaging method not recognized in DG_2_p1CG.\n");
    
    // Compute p1CG field at each DG node:
    for (Int ie = 0; ie < ne; ie++) {
        Int elemtype = mesh.elementype[ie];
        Int npv = mesh.elements[elemtype].npv;
        
        if (elemtype == 0) {
            double f1 = fieldVertices[mesh.t[0*ne+ie]];
            double f2 = fieldVertices[mesh.t[1*ne+ie]];
            double f3 = fieldVertices[mesh.t[2*ne+ie]];

            double a = -f1+f2;
            double b = -f1+f3;
            double d = f1;
            
            double f4, c;
            if (nd == 2)
                c = 0.0;
            else if (nd == 3) {
                f4 = fieldVertices[mesh.t[3*ne+ie]];
                c = -f1+f4;
            }

            for (Int j = 0; j < npv; j++) {
                double xi1 = ?? master.plocvl[0*npv+j];
                double xi2 = ?? master.plocvl[1*npv+j];
                double xi3;
                if (nd == 2)
                    xi3 = 0.0;
                else if (nd == 3)
                    ?? xi3 = master.plocvl[2*npv+j];

                if (activeElement[ie] == 1 || zeroPassiveElemets == 0)
                    p1CGfield[sol.scalarFieldStart[ie]+j] = a*xi1 + b*xi2 + c*xi3 + d;
                else
                    p1CGfield[sol.scalarFieldStart[ie]+j] = 0.0;
            }
        }
        else if (elemtype == 1) {
            double f1 = fieldVertices[mesh.t[0*ne+ie]];
            double f2 = fieldVertices[mesh.t[1*ne+ie]];
            double f3 = fieldVertices[mesh.t[2*ne+ie]];
            double f4 = fieldVertices[mesh.t[3*ne+ie]];
            double f5, f6, f7, f8;
            if (nd == 2) {
                f5 = 0.0;
                f6 = 0.0;
                f7 = 0.0;
                f8 = 0.0;
            }
            else if (nd == 3) {
                f5 = fieldVertices[mesh.t[4*ne+ie]];
                f6 = fieldVertices[mesh.t[5*ne+ie]];
                f7 = fieldVertices[mesh.t[6*ne+ie]];
                f8 = fieldVertices[mesh.t[7*ne+ie]];
            }

            for (Int j = 0; j < npv; j++) {
                double ?? xi1 = master.plocvl[0*npv+j];
                double ?? xi2 = master.plocvl[1*npv+j];
                double xi3;
                if (nd == 2)
                    xi3 = 0.0;
                else if (nd == 3)
                    ?? xi3 = master.plocvl[2*npv+j];
                
                double f12 = (1.0-xi1)*f1 + xi1*f2;
                double f43 = (1.0-xi1)*f4 + xi1*f3;
                double f1234 = (1.0-xi2)*f12 + xi2*f43;
                
                double f56 = (1.0-xi1)*f5 + xi1*f6;
                double f87 = (1.0-xi1)*f8 + xi1*f7;
                double f5678 = (1.0-xi2)*f56 + xi2*f87;

                if (activeElement[ie] == 1 || zeroPassiveElemets == 0)
                    p1CGfield[sol.scalarFieldStart[ie]+j] = (1.0-xi3)*f1234 + xi3*f5678;
                else
                    p1CGfield[sol.scalarFieldStart[ie]+j] = 0.0;
            }
        }
        else
            error("elemtype not recognized in DG_2_p1CG.\n");
    }
    
    delete[] numElementsAdj2Vertex; delete[] vertex2elem; delete[] vertex2node; delete[] avgFieldElem; delete[] fieldVertices; delete[] activeElement;
}

// // // // // void DG_2_p1CG(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims, Int avgMethod)
// // // // // {
// // // // //     Int nd = ndims[0];
// // // // // 
// // // // //     if (nd == 2)
// // // // //         DG_2_p1CG_2d(p1CGfield, DGfield, mesh, master, app, ndims, avgMethod);
// // // // //     else if (nd == 3)
// // // // //         DG_2_p1CG_3d(p1CGfield, DGfield, mesh, master, app, ndims, avgMethod);
// // // // //     else
// // // // //         error("Invalid number of dimensions.\n");
// // // // // }

#endif
