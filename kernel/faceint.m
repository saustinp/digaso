function [Ru, Rq, Rh, BD, F, E, G, H, Ju, Jq, Jh, wrb] = faceint(master,app,bf,dgnodes,UDG,UH,Ru,Rq,BD,Ju,Jq)
% FACEINTND compute face integrals 

% Obtain dgnodes, Jacobian matrix and determinant at Gauss points
[pg, nlg, jac] = facegeom(master.shapmf,dgnodes,master.permgeom);

ne   = size(UDG,2);
nd   = master.nd;
npv  = master.npv;
npf  = master.npf;
ngf  = master.ngf;
nfe  = size(master.perm,2);
nc   = app.nc;
ncu  = app.ncu;
nch  = app.nch;
arg  = app.arg;
time = app.time;
bcm  = app.bcm;        
bcs  = app.bcs;
bcd  = app.bcd;        
bcv  = app.bcv;
localsolve = app.localsolve;
adjoint = app.adjoint;
fbou   = str2func(app.fbou);
fhat   = str2func(app.fhat);

% Shap functions 
perm            = master.perm(:,:,1);
shapft          = master.shapft(:,:,1);
shapfg          = master.shapfg(:,:,1);
shapfgdotshapfc = master.shapfgdotshapfc(:,:,1);

if localsolve==0
    E = [];
else    
    tmp = zeros(npv,npf*nfe,ne);
    E = zeros(npv,npf*nfe,nd,ne);        
    for i=1:nd    
        njc = reshape(nlg(:,i).*(jac),[ngf,nfe*ne]);
        wrk = reshape(shapfgdotshapfc*njc,[npf npf*nfe ne]);   
        for j=1:nfe
            tmp(perm(:,j),(j-1)*npf+1:j*npf,:) = wrk(1:npf,(j-1)*npf+1:j*npf,:);      
        end    
        for k=1:ne
            E(:,:,i,k) = tmp(:,:,k);        
        end
    end
    for i = 1:nd
        Rq(:,:,i,:) = Rq(:,:,i,:) - reshape(mapContractK(E(:,:,i,:),UH,1,[2 3],4,1,3,2),[npv ncu 1 ne]);
    end
end
    
% DG solution at Gauss points
udgn = reshape(UDG(perm,:,:),[npf nfe*ne*nc]);
udgg = shapft*udgn;
udgg = reshape(udgg,[ngf*nfe*ne nc]);

% for comparison
udgn =reshape(udgn, [npf*nfe*ne nc]);

uh  = reshape(UH,[npf nfe*ne*nch]);
uhg = shapft*uh;
uhg = reshape(uhg,[ngf*nfe*ne nch]);

% Compute numerical fluxes 
[FH, FH_udg, FH_uh] = fhat(nlg, pg, udgg, uhg, arg, time);     

if isfield(app, 'debug_digaso')
    uhgdig=readbin('./debug/debug8.bin');  % uhg       
    ugfdig=readbin('./debug/debug9.bin');  % udgg      
    nlgfdig=readbin('./debug/debug10.bin'); % nlg      
    pgfdig=readbin('./debug/debug11.bin');  % pg       
    fhdig=readbin('./debug/debug12.bin');  % fh
    fhudig=readbin('./debug/debug13.bin');  % fh_u
    fhuhdig=readbin('./debug/debug14.bin');  % fh_uh
    ufdig = readbin('./debug/debug30.bin');
    udgdig=readbin('./debug/debug31.bin');  
            
    err = zeros(10,1);
    ugfdig = permute(reshape(ugfdig, [ngf nfe nc ne]),[1 2 4 3]);
    err(1) = max(abs(ugfdig(:)-udgg(:)));    
    
    ufdig = permute(reshape(ufdig, [npf nfe nc ne]),[1 2 4 3]);
    err(2) = max(abs(ufdig(:)-udgn(:)));
        
    udgdig = permute(reshape(udgdig, [npv nc ne]),[1 3 2]);        
    err(3) = max(abs(udgdig(:)-UDG(:)));    
    
    fprintf('|ugf-ugfdig| = %e,   |uf-ufdig| = %e,   |udg-udgdig| = %e  \n', err(1:3));     
    
    uhgdig = permute(reshape(uhgdig, [ngf nfe nch ne]),[1 2 4 3]);
    err(4) = max(abs(uhgdig(:)-uhg(:)));    
    
    pgfdig = permute(reshape(pgfdig, [ngf nfe nd ne]),[1 2 4 3]);
    err(5) = max(abs(pgfdig(:)-pg(:)));    
    
    nlgfdig = permute(reshape(nlgfdig, [ngf nfe nd ne]),[1 2 4 3]);
    err(6) = max(abs(nlgfdig(:)-nlg(:)));    
    
    fprintf('|uhg-uhgdig| = %e,   |pgf-pgfdig| = %e,   |nlg-nlgdig| = %e  \n', err(4:6));         
    
    fhdig = permute(reshape(fhdig, [ngf nfe nch ne]),[1 2 4 3]);
    tm = fhdig(:)-FH(:);
    err(7) = max(abs(tm));        
    
    fhudig = permute(reshape(fhudig, [ngf nfe nch nc ne]),[1 2 5 3 4]);
    tm = fhudig(:)-FH_udg(:);
    err(8) = max(abs(tm));        
    
    fhuhdig = permute(reshape(fhuhdig, [ngf nfe nch nch ne]),[1 2 5 3 4]);
    tm = fhuhdig(:)-FH_uh(:);
    err(9) = max(abs(tm));        
    
    fprintf('|fh-fhdig| = %e,   |fhu-fhudig| = %e,   |fhuh-fhuhdig| = %e  \n', err(7:9));         
    
    %error('Done comparing');  
end


% squeeze(pg)
% squeeze(nlg)
% squeeze(udgg)
% squeeze(uhg)
% squeeze(FH)
% pause
% tm = reshape(FH, [ngf nfe ne nch]);
% tm = reshape(tm(:,:,1,:),[ngf*nfe nch]);
% tn = reshape(FH_udg, [ngf nfe ne nch nc]);
% tn = reshape(tn(:,:,1,:,:),[ngf*nfe nch nc]);
% tp = reshape(FH_uh, [ngf nfe ne nch nch]);
% tp = reshape(tp(:,:,1,:,:),[ngf*nfe nch nch]);
% tm(:)
% squeeze(tn)
% [squeeze(tp) squeeze(pg(1:ngf*nfe,:))]
% arg{end}
% pause

% compute boundary fluxes
% fbou   = str2func(app.fbou);
% tm     = find(bf<0);
% nbf    = length(tm(:));
% JH    = zeros(ngf*nbf,ncu);
% JUDG   = zeros(ngf*nbf,nc);
% BH     = zeros(ngf*nbf,nch);
% BH_udg = zeros(ngf*nbf,nch,nc);
% BH_uh  = zeros(ngf*nbf,nch,nch);
% fbou   = str2func(app.fbou);
% an     = zeros(nbf*npf,1);
% bn     = zeros(nbf*ngf,1);
% j      = 1; 
% for k=1:ne                                
%     for is=1:nfe                            
%         if bf(is,k) < 0 % for each boundary face                        
%             im = (j-1)*ngf+1:j*ngf;            
%             in = (k-1)*nfe*ngf+(is-1)*ngf+1:(k-1)*nfe*ngf+is*ngf;            
%             bn(im) = in; 
%             ib = bcm(-bf(is,k));        
%             bv = repmat(bcs(-bf(is,k),:),[ngf 1]);            
%             [fh, dfh_dudg, dfh_duh] = fbou(ib,bv,nlg(in,:),pg(in,:),udgg(in,:),uhg(in,:),arg,time);                                                 
%             ia = bcd(-bf(is,k));        
%             bu = repmat(bcv(-bf(is,k),:),[ngf 1]);            
%             [~,djh_dudg,djh_duh] = fbou_adjoint(ia,bu,nlg(in,:),pg(in,:),udgg(in,:),uhg(in,:),arg,time);                           
%             JUDG(im,:) = djh_dudg;
%             JH(im,:) = djh_duh;
%             BH(im,:) = fh;
%             BH_udg(im,:,:) = dfh_dudg;
%             BH_uh(im,:,:) = dfh_duh;              
%             im = (j-1)*npf+1:j*npf;            
%             in = (k-1)*nfe*npf+(is-1)*npf+1:(k-1)*nfe*npf+is*npf;            
%             an(im) = in; 
%             j = j+1;
%         end                
%     end
% end

tm     = find(bf<0);
nbf    = length(tm(:));
BH     = zeros(ngf*nbf,nch);
BH_udg = zeros(ngf*nbf,nch,nc);
BH_uh  = zeros(ngf*nbf,nch,nch);
JH     = zeros(ngf*nbf,ncu);
JUDG   = zeros(ngf*nbf,nc);
an     = zeros(nbf*npf,1);
bn     = zeros(nbf*ngf,1);
j = 1;
for i = 1:length(bcm)    
    [I,J] = find(bf==-i);
    nfb = length(I);
    pgb  = zeros(ngf*nfb, size(pg,2));
    nlgb = zeros(ngf*nfb, nd);
    udgb = zeros(ngf*nfb, nc);
    uhgb = zeros(ngf*nfb, ncu);
    im = zeros(ngf*nfb,1);
    for k = 1:nfb
        i0 = (j-1)*ngf+1:j*ngf;
        i1 = (k-1)*ngf+1:k*ngf; 
        i2 = (J(k)-1)*nfe*ngf+(I(k)-1)*ngf+1:(J(k)-1)*nfe*ngf+I(k)*ngf; 
        bn(i0) = i2;
        im(i1)  = i0;
        pgb(i1,:)  = pg(i2,:);
        nlgb(i1,:) = nlg(i2,:);
        udgb(i1,:) = udgg(i2,:);
        uhgb(i1,:) = uhg(i2,:);
        i0 = (j-1)*npf+1:j*npf;
        i2 = (J(k)-1)*nfe*npf+(I(k)-1)*npf+1:(J(k)-1)*nfe*npf+I(k)*npf; 
        an(i0) = i2;
        j = j + 1;
    end          
    
    ib = bcm(i);
    bv = repmat(bcs(i,:),[ngf*nfb 1]);    
    if isempty(im)==0
    [fh, dfh_dudg, dfh_duh] = fbou(ib,bv,nlgb,pgb,udgb,uhgb,arg,time);                                                 
    BH(im,:) = fh;
    BH_udg(im,:,:) = dfh_dudg;
    BH_uh(im,:,:) = dfh_duh;              
    end
    
    if adjoint==1
        ia = bcd(i);        
        bu = repmat(bcv(i,:),[ngf*nfb 1]);            
        [~,djh_dudg,djh_duh] = fbou_adjoint(ia,bu,nlgb,pgb,udgb,uhgb,arg,time);                                               
        JUDG(im,:) = djh_dudg;
        JH(im,:) = djh_duh;
    end
end

if isfield(app, 'debug_digaso')
    fbdig=readbin('./debug/debug15.bin');  % fb
    fbudig=readbin('./debug/debug16.bin');  % fb_u
    fbuhdig=readbin('./debug/debug17.bin');  % fb_uh
    fbdig = reshape(fbdig,[ngf nch nbf]);
    fbudig = reshape(fbudig,[ngf nch nc nbf]);
    fbuhdig = reshape(fbuhdig,[ngf nch nch nbf]);
    
    FB     = zeros(ngf,nch,nbf);
    FB_udg = zeros(ngf,nch,nc,nbf);
    FB_uh  = zeros(ngf,nch,nch,nbf);    
    pg = reshape(pg, [ngf nfe ne nd]);
    nlg = reshape(nlg, [ngf nfe ne nd]);
    udgg = reshape(udgg, [ngf nfe ne nc]);
    uhg = reshape(uhg, [ngf nfe ne ncu]);
    im = 0;
    for ie = 1:size(bf,2)
      for j = 1:size(bf,1)
        if bf(j,ie)<0
          i = -bf(j,ie);
          ib = bcm(i);
          bv = repmat(bcs(i,:),[ngf 1]);
          pgb  = reshape(pg(:,j,ie,:),[ngf nd]);
          nlgb = reshape(nlg(:,j,ie,:),[ngf nd]);
          udgb = reshape(udgg(:,j,ie,:),[ngf nc]);
          uhgb = reshape(uhg(:,j,ie,:),[ngf ncu]);          
          
          if isempty(im)==0
          [fh, dfh_dudg, dfh_duh] = fbou(ib,bv,nlgb,pgb,udgb,uhgb,arg,time);                                                 
          im = im + 1;
          FB(:,:,im) = fh;
          FB_udg(:,:,:,im) = dfh_dudg;
          FB_uh(:,:,:,im) = dfh_duh;              
          e = fh - fbdig(:,:,im);
          if max(abs(e(:)))>1e-6
            [ie j i ib]
%              max(abs(e(:)))
%             fh
%             fbdig(:,:,im)
%             bcs(i,:)
%             pause
          end
          e = dfh_dudg - fbudig(:,:,:,im);
          if max(abs(e(:)))>1e-6
            [ie j i ib]
          end
          e = dfh_duh - fbuhdig(:,:,:,im);
          if max(abs(e(:)))>1e-6
            [ie j i ib]
          end
          
          end
        end
      end
    end
    
    err = zeros(10,1);
    tm = fbdig(:)-FB(:);
    err(7) = max(abs(tm));        
        
    tm = fbudig(:)-FB_udg(:);
    err(8) = max(abs(tm));        
    
    tm = fbuhdig(:)-FB_uh(:);
    err(9) = max(abs(tm));        
    
    fprintf('|fb-fbdig| = %e,   |fbu-fbudig| = %e,   |fbuh-fbuhdig| = %e  \n', err(7:9));         
    
    error('Done comparing');  
end

if adjoint==1
    wrk = reshape(bsxfun(@times,JH,jac(bn)),[ngf nbf*nch]);
    Jh = zeros(npf*nfe*ne,nch);
    Jh(an,:) = reshape(shapfg*wrk,[npf*nbf nch]);
    Jh = permute(reshape(Jh,[npf*nfe ne nch]), [3 1 2]);

    wrk = reshape(bsxfun(@times,JUDG,jac(bn)),[ngf nbf*nc]);
    Jun = zeros(npf*nfe*ne,nc);
    Jun(an,:) = reshape(shapfg*wrk,[npf*nbf nc]);
    Jun = reshape(Jun,[npf nfe ne nc]);
    Juq = zeros(npv, 1, ne, nc);
    for is = 1:nfe
       IP = perm(:,is);
       Juq(IP,1,:,:) = Juq(IP,1,:,:) + Jun(:,is,:,:); 
    end
    Juq = reshape(Juq,[npv ne nc]);
    Ju  = permute(Ju + Juq(:,:,1:ncu),[1 3 2]);
    Jq  = permute(Jq + Juq(:,:,ncu+1:end),[1 3 2]);
    %if nc>ncu
    %    Jq  = reshape(Jq, [npv ncu nd ne]);
    %end
else
    Ju = [];
    Jq = [];
    Jh = [];    
end

wrk = reshape(bsxfun(@times,FH,jac),[ngf nfe*ne*nch]);
Run = reshape(shapfg*wrk,[npf nfe ne nch]);

% jac
% reshape(FH,[ngf*nfe ne*nch])
% reshape(bsxfun(@times,FH,jac),[ngf*nfe ne*nch])
% reshape(Run,[npf*nfe ne*nch])
% pause

wrk = reshape(bsxfun(@times,FH_udg,jac),[ngf nfe*ne*nch*nc]);
BDn = reshape(shapfgdotshapfc*wrk,[npf npf nfe ne nch nc]);
wrb = wrk;

wrk = reshape(bsxfun(@times,FH_uh,jac),[ngf nfe*ne*nch*nch]);
Fn = reshape(shapfgdotshapfc*wrk,[npf npf nfe ne nch nch]);
   
Rut = zeros(npv,1,ne,ncu);
BDt = zeros(npv,npv,1,ne,ncu,nc);
F = zeros(npv,npf,nfe,ne,ncu,nch); 
for is=1:nfe  % Adding face contributions - vector dependencies avoided
    IP = perm(:,is);
    Rut(IP,1,:,:) = Rut(IP,1,:,:) + Run(:,is,:,:);
    BDt(IP,IP,1,:,:,:) = BDt(IP,IP,1,:,:,:) + BDn(:,:,is,:,:,:);
    F(IP,:,is,:,:,:) = F(IP,:,is,:,:,:) + Fn(:,:,is,:,:,:);
end
Ru = permute(Ru-reshape(Rut,[npv ne ncu]),[1 3 2]);
BD = permute(BD+reshape(BDt,[npv npv ne ncu nc]),[1 4 2 5 3]);
F = permute(reshape(F,[npv npf*nfe ne ncu nch]),[1 4 5 2 3]);

wrk = reshape(bsxfun(@times,BH,jac(bn)),[ngf nbf*nch]);
Run = reshape(Run,[npf*nfe*ne nch]);
Run(an,:) = reshape(shapfg*wrk,[npf*nbf nch]);
Rh  = -permute(reshape(Run,[npf*nfe ne nch]), [3 1 2]);

wrk = reshape(bsxfun(@times,BH_udg,jac(bn)),[ngf nbf*nch*nc]);
Gt = reshape(BDn,[npf npf*nfe*ne nch nc]);
Gt(:,an,:,:) = reshape(shapfgdotshapfc*wrk,[npf npf*nbf nch nc]);
Gt = reshape(Gt,[npf npf nfe ne nch nc]);

wrk = reshape(bsxfun(@times,BH_uh,jac(bn)),[ngf nbf*nch*nch]);
Ht = reshape(Fn,[npf npf*nfe*ne nch nch]);
Ht(:,an,:,:) = reshape(shapfgdotshapfc*wrk,[npf npf*nbf nch nch]);
Ht  = reshape(Ht,[npf npf nfe ne nch nch]);

G   = zeros(npv,npf,nfe,ne,nch,nc); 
H   = zeros(npf*nfe,npf*nfe,1,ne,nch,nch);
for is=1:nfe  
    IP = (is-1)*npf+1:is*npf;
    G(perm(:,is),:,is,:,:,:) = G(perm(:,is),:,is,:,:,:) + Gt(:,:,is,:,:,:);
    H(IP,IP,1,:,:,:) = H(IP,IP,1,:,:,:) + Ht(:,:,is,:,:,:);
end
G  = permute(reshape(G,[npv npf*nfe ne nch nc]),[4 2 1 5 3]);
H  = permute(reshape(H,[npf*nfe npf*nfe ne nch nch]),[4 1 5 2 3]);
% t  = reshape(H(:,:,:,:,1),[nch*npf*nfe nch*npf*nfe])
% pause



