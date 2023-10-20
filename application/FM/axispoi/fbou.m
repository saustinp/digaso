function [fh,fh_udg,fh_uh] = fbou(ib,ui,nl,p,udg,uh,param,time)
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
nch = 1;
nq = nc-nch;

tau   = param{end};
u     = udg(:,nch);
r = p(:,1);
% tau = r.*tau;   % Changed for axisymmetric

switch ib
    case 1  % Homogeneous dirichlet
        fh = tau.*(0-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = -tau;
    case 2  % Homogeneous neumann    
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        fh = fh + ui;
    case 3 % Inhomogeneous neumann
        l_ref = param{2};
        E_ref = param{4};
        phi0 = param{6};
        phi0_tilde = phi0/(E_ref*l_ref);

        fh = tau.*(phi0_tilde-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = -tau;

        % Below is the code for the axis poisson test case
%     case 1  % Dirichlet
%         % fh = tau.*(ui-uh);
%         % fh_udg = zeros(ng,nch,nc);
%         % fh_uh = -tau;

%         r = p(:,1);
%         y = p(:,2);
%         ui = exp(-y).*cos(r);

%         fh = tau.*(ui-uh);
%         fh_udg = zeros(ng,nch,nc);
%         fh_uh = -tau;
%     case 2  % "Extrapolate  m = u 
%         fh = tau.*(u-uh);
%         fh_u = tau;
%         fh_q = zeros(ng,nch,nq);
%         fh_udg = cat(3,fh_u,fh_q);
%         fh_uh = -tau;
%     case 3  % Prescribed flux
% %         x = p(:,1);
% %         y = p(:,2);
% %         ui = sin(x)*sin(y);        
%         [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
%         fh = fh + ui;
%     case 4  % Prescribed flux
%         x = p(:,1);        
%         %[p exp(-(x/0.4236).^4)]
%         qn = ui.*exp(-(x/0.4236).^4);
%         [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
%         fh = fh + qn;           
    otherwise
        error('unknown boundary type');
end

