function [fh,fh_udg,fh_uh] = fbou_poisson(ib,ui,nl,p,udg,uh,param,time)
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

r = p(:,1);
tau   = param{end};
tau = r.*tau;

u     = udg(:,nch);

switch ib
    case 1  % Dirichlet
        fh = tau.*(ui-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = -tau;
    case 2  % "Extrapolate  m = u 
        fh = tau.*(u-uh);
        fh_u = tau*ones(ng,nch,nch);
        fh_q = zeros(ng,nch,nq);
        fh_udg = cat(3,fh_u,fh_q);
        fh_uh = -tau*ones(ng,nch,nch);
    case 3  % Prescribed flux
%         x = p(:,1);
%         y = p(:,2);
%         ui = sin(x)*sin(y);        
        [fh,fh_udg,fh_uh] = fhat_poisson(nl,p,udg,uh,param,time);
        fh = fh + ui;
    case 4  % Prescribed flux
        x = p(:,1);        
        %[p exp(-(x/0.4236).^4)]
        qn = ui.*exp(-(x/0.4236).^4);
        [fh,fh_udg,fh_uh] = fhat_poisson(nl,p,udg,uh,param,time);
        fh = fh + qn; 
    case 5 % Dirichlet
        x = p(:,1);
        y = p(:,2);
        ui = sin(0.5*pi*x).*sin(0.5*pi*y);       
        fh = tau.*(ui-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = -tau;     
    otherwise
        error('unknown boundary type');
end


