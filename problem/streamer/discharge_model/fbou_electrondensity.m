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
% r = p(:,1);
r = 1;
tau = param{end};

switch ib
    case 1 % Electrode: species set to homogeneous nuemann and potential set to dirichlet
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);
        
        % Species densities: homogeneous neumann
        [fh_tmp,fh_udg_tmp,fh_uh_tmp] = fhat_axisymmetric(nl,p,udg,uh,param,time);

        fh(:,[1,2]) = fh_tmp(:,[1,2]);
        fh_udg(:,[1,2],:) = fh_udg_tmp(:,[1,2],:);
        fh_uh(:,[1,2],:) = fh_uh_tmp(:,[1,2],:);

        % Potential: dirichlet
        fh(:,3) = r.*tau .*(ui-uh(:,3));
        fh_uh(:,3,3) = -r.*tau;

    case 2  % Right "farfield"
        % Species + potential all have homogeneous neumann
        [fh,fh_udg,fh_uh] = fhat_axisymmetric(nl,p,udg,uh,param,time);

        % fh = zeros(ng,nch);
        % fh_udg = zeros(ng,nch,nc);
        % fh_uh = zeros(ng,nch,nch);

        % [fh_tmp,fh_udg_tmp,fh_uh_tmp] = fhat_axisymmetric(nl,p,udg,uh,param,time);

        % fh(:,1) = r.*tau .*(0-uh(:,1));
        % fh_uh(:,1,1) = -r.*tau;

        % fh(:,2) = r.*tau .*(0-uh(:,2));
        % fh_uh(:,2,2) = -r.*tau;

        % fh(:,3) = fh_tmp(:,3);
        % fh_udg(:,3,:) = fh_udg_tmp(:,3,:);
        % fh_uh(:,3,:) = fh_uh_tmp(:,3,:);

    case 3  % Symmetry boundary: extrapolate m=u or u=uhat
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);

        fh(:,1) = r.*tau.*(udg(:,1)-uh(:,1));
        fh_udg(:,1,1) = r.*tau;
        fh_uh(:,1,1) = -r.*tau;

        fh(:,2) = r.*tau.*(udg(:,2)-uh(:,2));
        fh_udg(:,2,2) = r.*tau;
        fh_uh(:,2,2) = -r.*tau;

        fh(:,3) = r.*tau.*(udg(:,3)-uh(:,3));
        fh_udg(:,3,3) = r.*tau;
        fh_uh(:,3,3) = -r.*tau;

    otherwise
        error('unknown boundary type');
end
