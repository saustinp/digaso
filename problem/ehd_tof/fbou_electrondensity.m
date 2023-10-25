function [fh,fh_udg,fh_uh] = fbou_electrondensity(ib,ui,nl,p,udg,uh,param,time)
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
%      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q

[ng,nc] = size(udg);
nch = nc/3;     % 3 components to UDG: (U,QX,QY) -> assuming 2D. For example, in 2D if U has 4 equations, nc will be 12

%%%%%%%%%       Copied from fluxgen_electrondensity.m
% Read in values from the p vector -> this normally only contains the (r,z) coordinates but we can also pass in additional fields as appended elements
r = p(:,1);

% Read in values from the u vector
% ne = udg(:,1);
% dne_dr = udg(:,2); % q is -grad(u)
% dne_dz = udg(:,3);
tau = param{end};

switch ib
    case 1  % Axisymmetry boundary
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);

        fh(:,1) = tau*(udg(:,1)-uh(:,1));
        fh_udg(:,1,1) = tau;
        fh_uh(:,1,1) = -tau;

    case 2  % Electrodes - homogeneous neumann
        [fh,fh_udg,fh_uh] = fhat_electrondensity(nl,p,udg,uh,param,time);

    case 3  % Farfield: homogeneous dirichlet
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);
        
        fh(:,1) = r.*tau .*(0-uh(:,1));
        fh_uh(:,1,1) = -r.*tau ;
    otherwise
        error('unknown boundary type');
end
