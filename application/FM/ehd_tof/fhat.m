% function [fh,fh_udg,fh_uh] = fhat_electrondensity(nl,p,udg,uh,param,time,signe)
% %FHAT flux function
% %   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
% %
% %      NL(N,ND)              Normal N points
% %      P(N,ND)               Coordinates for N points
% %      U(N,NC)               Unknown vector for N points with NC components
% %      Q(N,NC,ND)            Flux vector for N points with NC components in the
% %                            coordinate directions
% %      M(N,NC)               Hybrid unkowns for N points with NC components
% %      PARAM                 Parameter list
% %      FH(N,NC):              Volume flux at N points
% %      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
% %      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
% %      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q

% [ng,nc] = size(udg);
% nch = nc/3;     % 3 components to UDG: (U,QX,QY) -> assuming 2D. For example, in 2D if U has 4 equations, nc will be 12
% % nq = nc-nch;
% % nd = nq;

% % Physics parameters
% % r0 = param{1};
% % z0 = param{2};
% % s0 = param{3};
% % Nmax = param{4};
% % e = param{5};
% % epsilon0 = param{6};
% % Ua = param{7};
% % gamma = param{8};
% E_bd = param{9};
% r_tip = param{10};
% % n_ref = param{11};
% % N = param{12};
% mue_ref = param{13};
% tau = param{end};

% % These swarm parameters will eventually be replaced with calls to the curve fit functions.
% mu_e = .0378;
% mu_p = 2.43e-4;          % mu[3] Pos ion mobility [m^2/(Vs)]
% mu_n = 2.7e-4;           % mu[4] Neg mobility [m^2/(Vs)]
% De = 0.18;               % Electron diffusion coefficient [m^2/s]
% Dn = 0.043e-4;           % Neg diffusion coefficient [m^2/s]
% Dp = 0.028e-4;           % Pos ion diffusion coefficient [m^2/s]

% De_star = De/(mue_ref*E_bd*r_tip);
% Dn_star = Dn/(mue_ref*E_bd*r_tip);
% Dp_star = Dp/(mue_ref*E_bd*r_tip);

% % Note the u vector is: (u1, u2, u3, u4, q1x, q2x, q3x, q4x, q1y, q2y, q3y, q4y)
% %                        1   2   3    4     5       6       7     8     9       10      11    12
% % Note the u vector is: (ne, np, nn, phi, dne_dr, dnn_dr, dnp_dr, Er, dne_dz, dnn_dz, dnp_dz, Ez)

% % Without the electrostatic equation:
% % Note the u vector is: (u1, u2, u3, q1x, q2x, q3x, q1y, q2y, q3y)
% %                        1   2   3     4       5       6       7       8       9
% % Note the u vector is: (ne, np, nn, dne_dr, dnn_dr, dnp_dr, dne_dz, dnn_dz, dnp_dz)

% r = p(:,1);
% Ex = p(:,3);    % Ex is -grad(phi)
% Ey = p(:,4);

% % Read in values from the u vector
% ne = udg(:,1);
% nn = udg(:,2);
% np = udg(:,3);
% phi = udg(:,4);
% dne_dr = udg(:,5); % q is -grad(ne)
% dnn_dr = udg(:,6); % q is -grad(ne)
% dnp_dr = udg(:,7);
% Ex = udg(:,8);
% dne_dz = udg(:,9);
% dnn_dz = udg(:,10);
% dnp_dz = udg(:,11);
% Ey = udg(:,12);

% % Compute convective velocities for r and z components for each equation
% cr_e = -(mu_e/mue_ref)*Ex;
% cz_e = -(mu_e/mue_ref)*Ey;
% cr_n = -(mu_n/mue_ref)*Ex;
% cz_n = -(mu_n/mue_ref)*Ey;
% cr_p = (mu_p/mue_ref)*Ex;
% cz_p = (mu_p/mue_ref)*Ey;

% % ng x nch, same size as uh
% % This is the jacobian of the normal flux f_hat w.r.t UDG. UDG is [u, uq1x, q2x], and does _not_ involve uh.
% fh = zeros(ng,nch);
% fh(:,1) = r.*((cr_e.*uh(:,1)+De_star*dne_dr).*nl(:,1) + (cz_e.*uh(:,1)+De_star*dne_dz).*nl(:,2)) + tau.*(ne-uh(:,1));
% fh(:,2) = r.*((cr_n.*uh(:,2)+Dn_star*dnn_dr).*nl(:,1) + (cz_n.*uh(:,2)+Dn_star*dnn_dz).*nl(:,2)) + tau.*(nn-uh(:,2));
% fh(:,3) = r.*((cr_p.*uh(:,3)+Dp_star*dnp_dr).*nl(:,1) + (cz_p.*uh(:,3)+Dp_star*dnp_dz).*nl(:,2)) + tau.*(np-uh(:,3));

% % Jacobian of f_hat with respect to UDG.
% fh_udg = zeros(ng,nch,nc);
% fh_udg(:,1,1) = tau;                    % dfh(ne)_d(ne)
% fh_udg(:,1,4) = De_star*r.*nl(:,1);     % dfh(ne)_d(dne_dr)
% fh_udg(:,1,7) = De_star*r.*nl(:,2);     % dfh(ne)_d(dne_dz)

% fh_udg(:,2,2) = tau;                    % dfh(nn)_d(nn)
% fh_udg(:,2,5) = Dn_star*r.*nl(:,1);     % dfh(nn)_d(dnn_dr)
% fh_udg(:,2,8) = Dn_star*r.*nl(:,2);     % dfh(nn)_d(dnn_dz)

% fh_udg(:,3,3) = tau;                    % dfh(np)_d(np)
% fh_udg(:,3,6) = Dp_star*r.*nl(:,1);     % dfh(np)_d(dnp_dr)
% fh_udg(:,3,9) = Dp_star*r.*nl(:,2);     % dfh(np)_d(dnp_dz)

% % This is the jacobian of the normal flux f_hat w.r.t UHAT and does _not_ involve UDG.
% % For this problem the jacobian matrix is diagonal - each equation doesn't depend on uhat from the other equations
% fh_uh = zeros(ng,nch,nch);
% fh_uh(:,1,1) = cr_e.*r.*nl(:,1) + cz_e.*r.*nl(:,2) - tau;
% fh_uh(:,2,2) = cr_n.*r.*nl(:,1) + cz_n.*r.*nl(:,2) - tau;
% fh_uh(:,3,3) = cr_p.*r.*nl(:,1) + cz_p.*r.*nl(:,2) - tau;

% Copied from /HDG/application/ns/fhat.m
function [fh,fh_udg,fh_uh] = fhat_electrondensity(nl,p,udg,uh,param,time)
    %FHAT flux function
    %   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
    %
    %      NL(N,ND)              Normal N points
    %      P(N,ND)               Coordinates for N points
    %      U(N,NC)               Unknown vector for N points with NC components
    %      Q(N,NC,ND)            Flux vector for N points with NC components in the
    %                            coordinate directions
    %      M(N,NC)               Hybrid unkowns for N points with NC components
    %      PARAM                 Parameter list
    %      FH(N,NC):              Volume flux at N points
    %      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
    %      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
    %      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Lambda
    
    [ng,nch] = size(uh);
    
    % if nch==5
    %     [fh,fh_udg,fh_uh] = fhat3d(nl,p,udg,uh,param,time);
    %     return;
    % end
    
    nd = 2;    
        
    u = udg(:,1:nch);
    q = reshape(udg(:,nch+1:end),ng,nch,[]);
    
    % [f,f_udg] = flux(p,udg,param,time);
    % 
    % %[An,Anm] = getan(nl,uh,param,2);
    % An  = zeros(ng,nch,nch);
    % for i=1:nch
    %     An(:,i,i)=param{end};
    % end
    % Anm = zeros(ng,nch,nch,nch);
    % 
    % fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
    % %fh = permute(mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
    % 
    % fn_udg = mapContractK(f_udg,nl,[2 4],3,1,2,[],1);
    % 
    % fh_u = permute(fn_udg(:,1:nch,:),[3 1 2])+An;
    % fh_q = permute(fn_udg(:,nch+1:3*nch,:),[3 1 2]);
    % fh_udg = cat(3,fh_u,fh_q);
    
    % fh_uh = permute(mapContractK(Anm,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;   
    r = p(:,1); 
    
    [f,f_uhdg] = flux2d(p,[uh,udg(:,nch+1:end)],param,time);
    
    %[An,Anm] = getan(nl,uh,param,2);
    An  = zeros(ng,nch,nch);
    tau = param{end};
    for i=1:nch
        An(:,i,i)=r.*tau;
    end
    Anm = zeros(ng,nch,nch,nch);
    
    fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
    %fh = permute(mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
    
    fn_udgh = mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1);
    
    fh_u = An;
    fh_q = permute(fn_udgh(:,nch+1:(nd+1)*nch,:),[3 1 2]);
    fh_udg = cat(3,fh_u,fh_q);
    
    fh_uh = permute(fn_udgh(:,1:nch,:)+mapContractK(Anm,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;
    
    