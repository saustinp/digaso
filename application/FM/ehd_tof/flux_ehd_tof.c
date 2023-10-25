
#include "flux_ehd_tof2d.c"
#include "flux_ehd_tof3d.c"

void flux_ehd_tof(double *f, double *f_udg, double *pg, double *udg, double *odg, meshstruct &mesh,
        masterstruct &master, appstruct &app, double *param, double time, int ng, 
        int nc, int ncu, int nd, int ncd, int computeJacobian)
{     
    switch (nd) {
        case 2:            
            flux_ehd_tof2d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);            
            break;
        case 3:
            flux_ehd_tof3d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }    
}

