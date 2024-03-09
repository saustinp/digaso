function [Ru, Rq, BD, M, C, L, Q, Ju, Jq, wrl] = volint(master,app,dgnodes,UDG,SH,ODG)
    % VOLINTND compute volume integrals 
    
    % Obtain dgnodes, Jacobian matrix and determinant at Gauss points
    [pg, Xx, jac] = volgeom(master.shapmv,dgnodes);
    
    ne   = size(UDG,2);
    nd   = master.nd;
    npv  = master.npv;
    ngv  = master.ngv;
    
    nc   = app.nc;
    ncu  = app.ncu;
    arg  = app.arg;
    time = app.time;
    tdep = app.tdep;
    fc_q = app.fc_q;
    fc_u = app.fc_u;
    flg_p = app.flg_p;
    localsolve = app.localsolve;
    adjoint = app.adjoint;
    source = str2func(app.source);
    flux   = str2func(app.flux);
    
    % Shap functions and derivatives
    shapvt = master.shapvt;
    shapvg = reshape(master.shapvg,[npv ngv*(nd+1)]);
    shapvgdotshapvl = reshape(master.shapvgdotshapvl,[npv*npv ngv (nd+1)]);
    
    % DG solution at Gauss points
    udgg = reshape(UDG,[npv ne*nc]);
    udgg = shapvt(:,:,1)*udgg;
    udgg = reshape(udgg,[ngv*ne nc]);
    
    if isempty(ODG)==0
      nco   = size(ODG,3);
      odgg = reshape(ODG,[npv ne*nco]);
      odgg = shapvt(:,:,1)*odgg;
      odgg = reshape(odgg,[ngv*ne nco]);
    end
    
    if localsolve == 0    
        M  = [];
        C = [];
        L = [];
        Q = [];
        Rq = [];
    else
        % Mass matrix
        M = fc_q*reshape(shapvgdotshapvl(:,:,1)*reshape(jac,[ngv ne]),[npv npv ne]);
        
        % mass inverse times convection matrices
        C = zeros(npv,npv,nd,ne);
        for i=1:nd
            tmp = reshape(shapvgdotshapvl(:,:,2)*Xx(:,:,i,1),[npv npv ne]);
            for j=2:nd
                tmp = tmp + reshape(shapvgdotshapvl(:,:,j+1)*Xx(:,:,i,j),[npv npv ne]);
            end    
            for k=1:ne
                C(:,:,i,k) = tmp(:,:,k);        
            end
        end        
        
        Rq = reshape(mapContractK(M,SH(:,:,ncu+1:end)/fc_q-UDG(:,:,ncu+1:end),1,2,3,1,3,2),[npv ncu nd ne]);     
        for i=1:nd
            Rq(:,:,i,:) = Rq(:,:,i,:)+reshape(mapContractK(C(:,:,i,:),UDG(:,:,1:ncu),1,[2 3],4,1,3,2),[npv ncu 1 ne]);
        end
        
        if flg_p
            ncq = nc - ncu - 1;        
            [w, w_q] = press(udgg(:,ncu+2:end)); 
            L = shapvg(1:npv,1:ngv)*reshape(w.*jac,[ngv ne]); % npv x ne
            Q = shapvgdotshapvl(:,:,1)*reshape(bsxfun(@times,w_q,reshape(jac,[ngv*ne 1])),[ngv ne*ncq]);
            Q = permute(reshape(Q,[npv npv ne ncq]),[1 2 4 3]);
        else
            L = [];
            Q = [];
        end    
            
    %     Rq(:,:,1,1)
    %     Rq(:,:,2,1)
    %     pause
    end
    % udgg(1:ngv,:)
    % pause
    
    % d ni/ nt = s_e -> (fc_u * ni - odg) = s_e -> ni = (1/fc_u) * (odg + s_e) 
    % Laplacian (phi) = (ni - ne)
    % Laplacian (phi) = ((1/fc_u) * (odg + s_e(ne, E)) - ne)
    
    
    % Fluxes and source at Gauss points   
    if isempty(ODG)==0
      [f, f_udg] = flux( pg, udgg, odgg, arg, fc_u, time);
      [s, s_udg] = source( pg, udgg, odgg, arg, fc_u, time); 
    else
      [f, f_udg] = flux( pg, udgg, arg, time);
      [s, s_udg] = source( pg, udgg, arg, time);       
    end
    
    % Code to test the C->Matlab flux and source functions
    % Flag implemented so we don't have to uncomment it every time
    if isfield(app,'check_volint_flg')
    
        % Prepare data structures
        nch=ncu;
        ncd=3;
        [ng,nc] = size(udgg);
        arg2 = cell2mat(arg);
        pg_digaso = reshape(pg,1,[]);   % flatten
        udg_digaso = reshape(udgg,1,[]);    % flatten
        s_digaso = zeros(size(s));
        s_digaso = reshape(s_digaso,1,[]);
        s_udg_digaso = zeros(size(s_udg));
        s_udg_digaso = reshape(s_udg_digaso,1,[]);
        f_digaso = zeros(size(f));
        f_digaso = reshape(f_digaso,1,[]);
        f_udg_digaso = zeros(size(f_udg));
        f_udg_digaso = reshape(f_udg_digaso,1,[]);
    
        a1=zeros(1,size(s_digaso,2));
        a2=zeros(1,size(s_udg_digaso,2));
        a3=zeros(1,size(pg_digaso,2));
        a4=zeros(1,size(udg_digaso,2));
        a5=zeros(1,size(arg2,2));
    
        disp('Generating MEX file for source ...')
        codegen source_digaso_matlab -args {a1,a2,a3,a4,a5,1,1,1,1,1,1} ../digaso_compare/source.c
    
        a1=zeros(1,size(f_digaso,2));
        a2=zeros(1,size(f_udg_digaso,2));
        disp('Generating MEX file for flux ...')
        codegen flux_digaso_matlab -args {a1,a2,a3,a4,a5,1,1,1,1,1,1} ../digaso_compare/flux.c
    
        % Run source function check
        [s_digaso, s_udg_digaso] = source_digaso_matlab_mex(s_digaso, s_udg_digaso, pg_digaso, udg_digaso, arg2, time, ng, nc, ncu, nd, ncd);
        s_digaso = reshape(s_digaso, [ng,nch]);
        s_udg_digaso = reshape(s_udg_digaso, [ng,nch,nc]);
    
        % Run flux function check
        [f_digaso, f_udg_digaso] = flux_digaso_matlab_mex(f_digaso, f_udg_digaso, pg_digaso, udg_digaso, arg2, time, ng, nc, ncu, nd, ncd);
        f_digaso = reshape(f_digaso, [ng,nch,nd]);
        f_udg_digaso = reshape(f_udg_digaso, [ng,nch,nd,nc]);
        
        str=sprintf('Linf norm in s: %d', max(max(s-s_digaso)));
        disp(str)
        str=sprintf('Linf norm in s_udg: %d', max(max(max(s_udg-s_udg_digaso))));
        disp(str)
        str=sprintf('Linf norm in f: %d', max(max(max(f-f_digaso))));
        disp(str)
        str=sprintf('Linf norm in f_udg: %d', max(max(max(max(f_udg-f_udg_digaso)))));
        disp(str)
    
        error('Done checking volume integral ...')
    
    end
    
    % Check FLUX with the C++ code
    % disp('DIGASO')
    % f_dig = readbin('/Users/saustin/Dropbox (MIT)/DIGASOV2/problem/streamer/discharge_model/debug/debug4.bin');
    % f_udg_dig = readbin('/Users/saustin/Dropbox (MIT)/DIGASOV2/problem/streamer/discharge_model/debug/debug5.bin');
    % 
    % a=sort([f(:) f_dig],1);
    % max(abs(a(:,1)-a(:,2)))
    % a=sort([f_udg(:) f_udg_dig],1);
    % max(abs(a(:,1)-a(:,2)))
    
    f     = reshape(f,[ngv ne ncu nd]);
    f_udg = permute(reshape(f_udg,[ngv ne ncu nd nc]),[1 2 3 5 4]); 
    s     = reshape(s(:,1:ncu),[ngv*ne ncu]);
    s_udg = reshape(s_udg(:,1:ncu,1:nc),[ngv*ne ncu nc]);
    
    
    if adjoint==1
        [~,JUDG] = source_adjoint(pg, udgg, arg, time);     
        JUDG = bsxfun(@times,reshape(JUDG(:,:,:,1),[ngv ne nc]),reshape(jac,[ngv ne 1]));        
        Juq = reshape(master.shapvg(:,:,1)*reshape(JUDG,[ngv ne*nc]),[npv ne nc]);                 
        Ju = Juq(:,:,1:ncu);
        Jq = Juq(:,:,ncu+1:end);
    else
        Ju = [];
        Jq = [];
    end
    
    % Update source term for time-dependent problems
    if tdep    
        Stn = reshape(SH(:,:,1:ncu),[npv ne*ncu]);
        Stg = shapvt(:,:,1)*Stn;
        Stg = reshape(Stg,[ngv*ne ncu]);
    
        % axis symmetry
        axisflag = isfield(app,'axisymmetry');
        if axisflag
            xr = app.axisymmetry.*abs(pg(:,1));     % pg=the first component of DGnodes (x or r coord)
            xr = reshape(xr,[ngv*ne 1]);
        else
            xr = ones(ngv*ne,1);
        end
        
        % time-derivative coefficients
        dtcoef = isfield(app,'dtcoef');
        if dtcoef
            dtcoef = app.dtcoef;
        else
            dtcoef = ones(ncu,1);
        end

    %    s = s + bsxfun(@times,xr,Stg-udgg(:,1:ncu)*fc_u);
    %     s = s + Stg - udgg(:,1:ncu)*fc_u;    
        % axis symmetry
        %fc_u
        
        fcu_vector = isfield(app,'fcu_vector');
        if fcu_vector
            dtcoef = dtcoef.*app.fcu_vector;
            for i=1:ncu
                s(:,i) = s(:,i) + dtcoef(i)*xr.*(Stg(:,i)-udgg(:,i)*fc_u);
                s_udg(:,i,i) = s_udg(:,i,i) - dtcoef(i)*fc_u*xr;
            end    
        else
            for i=1:ncu
                s(:,i) = s(:,i) + dtcoef(i)*xr.*(Stg(:,i)-udgg(:,i)*fc_u);
                s_udg(:,i,i) = s_udg(:,i,i) - dtcoef(i)*fc_u*xr;
        %         s_udg(:,i,i) = s_udg(:,i,i) - fc_u;
            end    
        end
    end    
    
    if isfield(app, 'debug_digaso')
        fdig=readbin('./debug/debug4.bin');
        fudig=readbin('./debug/debug5.bin');
        sdig=readbin('./debug/debug6.bin');
        sudig=readbin('./debug/debug7.bin');

        fdig = permute(reshape(fdig,[ngv ncu nd ne]),[1 4 2 3]);
        fudig = permute(reshape(fudig,[ngv ncu nd nc ne]),[1 5 2 4 3]);
        sdig = permute(reshape(sdig,[ngv ncu ne]),[1 3 2]);
        sudig = permute(reshape(sudig,[ngv ncu nc ne]),[1 4 2 3]);
                
        err = zeros(10,1);
        err(1) = max(abs(f(:)-fdig(:)));
        err(2) = max(abs(f_udg(:)-fudig(:)));
        err(3) = max(abs(s(:)-sdig(:)));
        err(4) = max(abs(s_udg(:)-sudig(:)));
        fprintf('|f-fdig| = %e,   |fu-fudig| = %e  \n', err(1:2));
        fprintf('|s-sdig| = %e,   |su-sudig| = %e  \n', err(3:4));
    end
    
    % compute wrk and wrl to time with shape functions
    wrk = zeros(ngv*(nd+1),ne*ncu);
    wrl = zeros(ngv*(nd+1),ne*ncu*nc);
    wrk(1:ngv,:) =  reshape(bsxfun(@times,s,jac),[ngv ne*ncu]);
    wrl(1:ngv,:) = -reshape(bsxfun(@times,s_udg,reshape(jac,[ngv*ne 1 1])),[ngv ne*ncu*nc]);
    for i=1:nd
        fk = bsxfun(@times,f(:,:,:,1),Xx(:,:,1,i));
        fl = bsxfun(@times,f_udg(:,:,:,:,1),Xx(:,:,1,i));
        for j=2:nd
            fk = fk + bsxfun(@times,f(:,:,:,j),Xx(:,:,j,i));
            fl = fl + bsxfun(@times,f_udg(:,:,:,:,j),Xx(:,:,j,i));
        end
        wrk(i*ngv+1:(i+1)*ngv,:) = reshape(fk,[ngv ne*ncu]);
        wrl(i*ngv+1:(i+1)*ngv,:) = -reshape(fl,[ngv ne*ncu*nc]);
    end
    
    % Volume residual
    % [Phi Phi_xi Phi_eta] x [S.*jac; Fx.*Xx(:,:,1,1)+Fy.*Xx(:,:,2,1); Fx.*Xx(:,:,1,2)+Fy.*Xx(:,:,2,2)]
    Ru = shapvg*wrk; % [npv ngv*(nd+1)] x [ngv*(nd+1) ne*ncu] 
    Ru = reshape(Ru,[npv ne ncu]); 
    
    % volume matrices
    BD = reshape(shapvgdotshapvl,[npv*npv ngv*(nd+1)])*wrl;
    BD = reshape(BD,[npv npv ne ncu nc]);
    

    if isfield(app, 'debug_digaso')
        Mdig=readbin('./debug/debug1.bin');
        Cdig=readbin('./debug/debug2.bin');
        BDdig=readbin('./debug/debug18.bin');

        BD1 = permute(reshape(BD(:,:,1,:,:), [npv npv ncu nc]), [1 3 2 4]);
        B = reshape(BD1(:,:,:,ncu+1:end), [npv ncu npv ncu nd]);
        D = reshape(BD1(:,:,:,1:ncu), [npv ncu npv ncu]);
        MiC = zeros(npv,npv,nd);
        Mi = inv(M(:,:,1));
        for d = 1:nd
            MiC(:,:,d) = Mi*C(:,:,d,1);
        end
        BMiC = mapContractK(B, MiC,[1 2 4],[3 5],[],[1 3],2,[]);
        
        BMiC = permute(BMiC,[1,2,4,3]);
        D = D + BMiC;
        BDdig = BDdig(1:npv*ncu*npv*ncu);
        max(abs(D(:) - BDdig(:)));
        

        % BD_dig = permute(BD_dig,[1,3,4,2]);
        % npv npv ne ncu nc
        
        % M_matlab = M(:,:,1);
        % disp('MATLAB')
        % M_matlab = chol(M_matlab);
        % Rq_mat = Rq(:,:,:,1);
        % Ru_mat = Ru(:,1,:);
        % BD_mat = BD(:,:,1,:,:);
        % 
        % disp('Rq')
        % max(abs(Rq_dig(:) - Rq_mat(:)))
        % disp('Ru')
        % max(abs(Ru_dig(:) - Ru_mat(:)))
        % disp('BD')
        % max(abs(BD_dig(:) - BD_mat(:)))
        % disp('dgnodes')
        % dgnodes = dgnodes(:,1,:);
        % max(abs(dgnodes_dig(:) - dgnodes(:)))
        % 
        err = zeros(10,1);
        err(1) = max(abs(Mdig-M(:)));
        err(2) = max(abs(Cdig-C(:)));
        % err(2) = max(abs(dig-(:)));
        fprintf('|M-Mdig| = %e,   |C-Cdig| = %e,   |BD-BDdig| = %e  \n', err(1:3));
    end
    
    
    
    % % F: nge*ne*ncu*nd
    % % S: nge*ne*ncu
    % % Xx: nge*ne*nd*nd
    
    % jac = reshape(jac,[ngv ne]);
    % [f, f_udg] = flux( pg, udgg, arg, time);
    % [s, s_udg] = source( pg, udgg, arg, time); 
    % f     = reshape(f,[ngv ne ncu nd]);
    % f_udg = reshape(f_udg,[ngv ne ncu nd nc]);
    % s     = reshape(s(:,1:ncu),[ngv ne ncu]);
    % s_udg = reshape(s_udg(:,1:ncu,1:nc),[ngv ne ncu nc]);
    % 
    % wrk = zeros(ngv,(nd+1),ncu,ne);
    % for i = 1:ngv
    %     for j = 1:ne
    %         for k = 1:ncu
    %             wrk(i,1,k,j) = s(i,j,k)*jac(i,j);
    %             for m = 1:nd
    %                 wrk(i,m+1,k,j) = f(i,j,k,1)*Xx(i,j,1,m);
    %                 for n = 2:nd
    %                     wrk(i,m+1,k,j) = wrk(i,m+1,k,j) + f(i,j,k,n)*Xx(i,j,n,m);
    %                 end
    %             end
    %         end
    %     end
    % end
    % Rut = shapvg*reshape(wrk,[ngv*(nd+1) ncu*ne]); % [npv ngv*(nd+1)] x [ngv*(nd+1) ncu*ne] 
    % Rut = reshape(Rut,[npv ncu ne]); 
    % e = Ru - permute(Rut,[1 3 2]);
    % max(abs(e(:)))
    % 
    % %f_udg = reshape(f_udg,[ngv ne ncu nd nc]);
    % wrk = zeros(ngv,(nd+1),ncu,nc,ne);
    % for i = 1:ngv
    %     for j = 1:ne
    %         for k = 1:ncu
    %             for l = 1:nc
    %                 wrk(i,1,k,l,j) = s_udg(i,j,k,l)*jac(i,j);
    %                 for m = 1:nd
    %                     wrk(i,m+1,k,l,j) = f_udg(i,j,k,1,l)*Xx(i,j,1,m);
    %                     for n = 2:nd
    %                         wrk(i,m+1,k,l,j) = wrk(i,m+1,k,l,j) + f_udg(i,j,k,n,l)*Xx(i,j,n,m);
    %                     end
    %                 end
    %             end
    %         end
    %     end
    % end
    % BDt = reshape(-shapvgdotshapvl,[npv*npv ngv*(nd+1)])*reshape(wrk,[ngv*(nd+1) ncu*nc*ne]);
    % BDt = reshape(BDt,[npv npv ncu nc ne]); 
    % e = BD - permute(BDt,[1 2 5 3 4]);
    % max(abs(e(:)))
    % 
    % if localsolve ~= 0 
    %     Mt = fc_q*reshape(shapvgdotshapvl(:,:,1)*jac,[npv npv ne]);
    %     wrk = permute(Xx,[1 4 3 2]);
    %     Ct = reshape(shapvgdotshapvl(:,:,2:end),[npv*npv ngv*nd])*reshape(wrk,[ngv*nd nd*ne]);
    %     Ct = reshape(Ct,[npv npv nd ne]); 
    %     e = C - Ct;
    %     max(abs(e(:)))
    % end
    % 
    