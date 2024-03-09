% Copied from /HDG/application/ns/fhat.m
function [fh,fh_udg,fh_uh] = fhat2_axisymmetric(nl,p,udg,odg,uh,param,factor,time)
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
    
    nd = 2;    
        
    u = udg(:,1:nch);
    q = reshape(udg(:,nch+1:end),ng,nch,[]);
    
    r = p(:,1); 
    
    [f,f_uhdg] = flux2d2(p,[uh,udg(:,nch+1:end)],odg,param,factor,time);
    
    An  = zeros(ng,nch,nch);
    tau = param{end};
    for i=1:nch
        % An(:,i,i)=r.*tau;
        An(:,i,i)=tau;
    end
    Anm = zeros(ng,nch,nch,nch);
    
    fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
    
    fn_udgh = mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1);
    
    fh_u = An;
    fh_q = permute(fn_udgh(:,nch+1:(nd+1)*nch,:),[3 1 2]);
    fh_udg = cat(3,fh_u,fh_q);
    
    fh_uh = permute(fn_udgh(:,1:nch,:)+mapContractK(Anm,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;
    
    