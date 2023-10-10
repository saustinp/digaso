function [fh,fh_udg,fh_uh] = fboucorona(ib,ui,nl,p,udg,uh,param,time)
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
nch = size(uh,2);
nq = nc-nch;

kappa = param{4};
tau   = param{end};
u     = udg(:,1:nch);

switch ib
    case 1  % Dirichlet      
        tau=1;
        if (max(p(:,1))>8000)            
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = tau*(ui(:,1)-uh(:,1));
        fh(:,2) = tau*(u(:,2)-uh(:,2));   
        fh_udg(:,2,2) = tau;        
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end            
        else                        
        fh = tau.*(ui(:,1:nch)-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end            
        end
        
%         if (size(p,1)==0) || (max(p(:,1))<3500)
%         fh = tau.*(ui(:,1:nch)-uh);
%         fh_udg = zeros(ng,nch,nc);
%         fh_uh =  zeros(ng,nch,nch);
%         for i = 1:nch           
%             fh_uh(:,i,i) = -tau;
%         end
%         else            
%         fh =  zeros(ng,nch);
%         fh_udg = zeros(ng,nch,nc);
%         fh_uh =  zeros(ng,nch,nch);
%         fh(:,1) = tau*(ui(:,1)-uh(:,1));
%         fh(:,2) = tau*(u(:,2)-uh(:,2));   
%         fh_udg(:,2,2) = tau;        
%         for i = 1:nch           
%             fh_uh(:,i,i) = -tau;
%         end            
%         end
    case 2  % "Extrapolate  m = u 
        fh = tau*(u-uh);
        fh_u = zeros(ng,nch,nch);        
        fh_q = zeros(ng,nch,nq);
        fh_uh = zeros(ng,nch,nch);
        for i = 1:nch
            fh_u(:,i,i) = tau;
            fh_uh(:,i,i) = -tau;
        end
        fh_udg = cat(3,fh_u,fh_q);        
    case 3  % dirichlet and extrapolation        
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = tau*(ui(:,1)-uh(:,1));
        fh(:,2) = tau*(u(:,2)-uh(:,2));   
        fh_udg(:,2,2) = tau;        
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end
    case 4 % fluxes
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        fh = fh + ui;    
    case 5 % symmetry
        nx = nl(:,1);
        ny = nl(:,2);
        phi = udg(:,1);
        rho = udg(:,2);        
        phix = udg(:,3);
        rhox = udg(:,4);
        phiy = udg(:,5);
        rhoy = udg(:,6);
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = phix.*nx + phiy.*ny + tau*(phi - uh(:,1));
        fh(:,2) = rhox.*nx + rhoy.*ny + tau*(rho - uh(:,2));        
        fh_udg(:,1,1) = tau;
        fh_udg(:,1,3) = nx;
        fh_udg(:,1,5) = ny;
        fh_udg(:,2,2) = tau;
        fh_udg(:,2,4) = nx;
        fh_udg(:,2,6) = ny;        
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end
    case 6 % plate
        nx = nl(:,1);
        ny = nl(:,2);        
        rho = udg(:,2);                
        rhox = udg(:,4);        
        rhoy = udg(:,6);
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = tau*(ui(:,1) - uh(:,1));
        fh(:,2) = rhox.*nx + rhoy.*ny + tau*(rho - uh(:,2));        
        fh_udg(:,2,2) = tau;
        fh_udg(:,2,4) = nx;
        fh_udg(:,2,6) = ny;        
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end        
    case 7  % Dirichlet        
        fh = tau.*(ui(:,1:nch)-uh);
        load('boundata3.mat');
        xc = p(:,1);
        yc = p(:,2);
        t1 = cart2pol((xc-xw),(yc-yw));
        t1(t1<0) = t1(t1<0)+2*pi;
        b1 = interp1(etaw,Peek,t1);
        u1 = b1-R;
        u1(u1<0) = 0.0;
        u1 = kappa*u1;
        fh(:,2) = tau.*(u1-uh(:,2));        
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end        
    case 8  % Dirichlet             
        xw = param{5};
        yw = param{6};
        t1 = param{7};
        t2 = param{8};        
        fh = tau.*(ui(:,1:nch)-uh);
        xc = p(:,1);
        yc = p(:,2);
        tt = cart2pol((xc-xw),(yc-yw));
        tt(tt<0) = tt(tt<0)+2*pi;
        tm = 0.5*(t1+t2);        
        tc = (log(1e-3))/(abs(t1-tm).^3);        
        ft = exp(tc*abs(tt-tm).^3);
        ft(tt>t2)=0;
        ft(tt<t1)=0;
        u1 = kappa*ft;        
        fh(:,2) = tau.*(u1-uh(:,2));        
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end                
    otherwise
        error('unknown boundary type');
end

