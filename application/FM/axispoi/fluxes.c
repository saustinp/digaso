#include "flux_axispoisson2d.c"
#include "flux_axispoisson3d.c"

// Written by: C. Nguyen & P. Fernandez

void flux_axispoisson(double *f, double *f_udg, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            flux_axispoisson2d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            flux_axispoisson3d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

void fluxonly_axispoisson(double *f, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fluxonly_axispoisson2d(f, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fluxonly_axispoisson3d(f, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}


#include "source_axispoisson2d.c"
#include "source_axispoisson3d.c"

void source_axispoisson(double *s, double *s_udg, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            source_axispoisson2d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            source_axispoisson3d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

void sourceonly_axispoisson(double *s, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            sourceonly_axispoisson2d(s, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            sourceonly_axispoisson3d(s, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

#include "fhat_axispoisson2d.c"
#include "fhat_axispoisson3d.c"

void fhat_axispoisson(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, 
        double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fhat_axispoisson2d(fh, fh_udg, fh_uh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fhat_axispoisson3d(fh, fh_udg, fh_uh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        default:
            exit(-1);
            break;
    }
}

void fhatonly_axispoisson(double *fh, double *pg, double *udg, 
        double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fhatonly_axispoisson2d(fh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fhatonly_axispoisson3d(fh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        default:
            exit(-1);
            break;
    }
}

#include "fbou_axispoisson2d.c"
#include "fbou_axispoisson3d.c"

void fbou_axispoisson(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nch, int nc, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fbou_axispoisson2d(fh, fh_u, fh_uh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;        
        case 3:            
            fbou_axispoisson3d(fh, fh_u, fh_uh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}


void fbouonly_axispoisson(double *fh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nch, int nc, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fbouonly_axispoisson2d(fh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;        
        case 3:            
            fbouonly_axispoisson3d(fh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

