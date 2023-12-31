void fbou_ehd_tof2d(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nc, int nch, int nd, int ncd)
{            
    
    int    i;                   
    double kappa   = param[0];
    double tau = param[1];   
    double r;       // Mod for axisymmetric
    
    if (ib==1) { /* Axisymmetry boundary */
        for (i=0; i<ng*nch*nc; i++)
            fh_u[i] = 0.0;
        for (i=0; i<ng*nch*nch; i++)
            fh_uh[i] = 0.0;

        // Can change the way you loop through these arrays to have the desired effect
        for (i=0; i<ng*nch; i++) {
            fh[i] = tau*(udg[i]-uhg[i]);    // Only need the first component of udg. hence the udg[i].
            fh_u[i] = tau;                      // Only need the first component of fh_u, hence the fh_u[i]
        }

        // for (i=0; i<ng*nch*nc; i++) {   
        //     fh_u[i] = 0.0;
        // }
        for (i=0; i<ng*nch*nch; i++){
            fh_uh[i] = -tau;                // In this case, nch=1 so this is the same size as fh
        }          
    }                                               
    else if (ib==2) { /* Neumman */
        fhat_ehd_tof2d(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd);    
        // for (i=0; i<ng*nch; i++) 
        //     fh[i] = fh[i] + ui[0];
    }                                     
    else if (ib==3) { /* Homogeneous dirichlet */
        for (i=0; i<ng*nch; i++) {
            fh[i] = tau*(0-uhg[i]);
        }
        for (i=0; i<ng*nch*nc; i++) {   
            fh_u[i] = 0.0;
        }
        for (i=0; i<ng*nch*nch; i++){
            fh_uh[i] = -tau;                // In this case, nch=1 so this is the same size as fh
        }          
    }                           
    else {                        
        printf("This error is in %s on line %d\n",__FILE__, __LINE__);
        printf("Boundary condition %d is not implemented yet.", ib);            
        exit(-1);                                    
    }                
}


// void fbouonly_poisson2d(double *fh, double *pg, double *udg, double *uhg, double *nl,
//           double *ui, double *param, double time, int ib,
//           int ng, int nc, int nch, int nd, int ncd)
// {        
//     
//     int    i;                   
//     double kappa   = param[0];
//     double tau = param[1];    
//     
//     if (ib==1) { /* Dirichlet */
//         for (i=0; i<ng*nch; i++) 
//             fh[i] = tau*(ui[0]-uhg[i]);           
//     }                                     
//     else if (ib==2) { /* Neumman */
//         fhatonly_poisson(fh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
//         
//         for (i=0; i<ng*nch; i++) 
//             fh[i] = fh[i] + ui[0];
//     }                                     
//     else if (ib==3) { /* Neumman */
//         fhatonly_poisson(fh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
//         
//         for (i=0; i<ng*nch; i++) 
//             fh[i] = fh[i] + nl[i];
//     }                                     
//     else if (ib==4) { /* Neumman */
//         fhatonly_poisson(fh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
//         
//         for (i=0; i<ng*nch; i++) 
//             fh[i] = fh[i] + nl[ng+i];
//     }                                     
//     else if (ib==5) { /* Neumman */
//         fhatonly_poisson(fh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
//         
//         for (i=0; i<ng*nch; i++) 
//             fh[i] = fh[i] + nl[2*ng+i];
//     }                                     
//     else {                        
//         printf("This error is in %s on line %d\n",__FILE__, __LINE__);
//         printf("Boundary condition %d is not implemented yet.", ib);            
//         exit(-1);                                    
//     }                
// }
// 
// 
