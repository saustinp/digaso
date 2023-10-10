#ifndef __GMRES_H
#define __GMRES_H

#include "../utilities/wrapper.h"
#include "../utilities/datastructs.h"

void computeOrdering(sysstruct &sys, Int preconditioner);

void computePreconditioner(sysstruct &sys, Int preconditioner);

void matvec(double* v, sysstruct &sys, double* x);

void matvecMPI(double* v, sysstruct &sys, double* x, Int sendInputVector, Int sendOutputVector, Int * requestCounter, double* matvectimes);

void gmres(sysstruct &sys, Int *convFlag);

void gmresmpi(sysstruct &sys, Int *convFlag);

#endif
