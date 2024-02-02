% These inputs must match the dimensions of the input vectors in the final output function. For example:
% function [f,f_udg] = flux_electrondensity(p,udg,param,time)

syms x1 x2  % Variables stored in the p vector
syms u1 u2 u3 u4 u5 u6 u7 u8 u9     % Components of UDG
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10    % Components of the physics parameters
syms zero one

param = [param1 param2 param3 param4 param5 param6 param7 param8 param9 param10];
udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9];
p = [x1 x2];

% Read in values from the p vector
r_tilde = p(1);

% Read in values from the u vector
ne_tilde = udg(1);
ni_tilde = udg(2);
phi_tilde = udg(3);
dne_dr_tilde = udg(4);
dni_dr_tilde = udg(5);
Er_tilde = udg(6);
dne_dz_tilde = udg(7);
dni_dz_tilde = udg(8);
Ez_tilde = udg(9);

% Load physics param
l_ref = param(1);
mu_ref = param(2);
E_ref = param(3);

% Compute transport coefficients
normE_tilde = sqrt(Er_tilde^2 + Ez_tilde^2);
normE = normE_tilde*E_ref;

mue_tilde = (2.3987*normE^(-.26))/mu_ref;
De_tilde = 4.3628e-3*normE^.22 / (l_ref*mu_ref*E_ref);

fv = [De_tilde.*dne_dr_tilde,0,Er_tilde,...
      De_tilde.*dne_dz_tilde,0,Ez_tilde];       % The negative sign is included in eqns 1 and 3 bc q=-grad(u)

fi = [-Er_tilde*mue_tilde*ne_tilde,0,0,...
      -Ez_tilde*mue_tilde*ne_tilde,0,0];

% Multiply flux by r for axisymmetry
f1 = r_tilde*(fi + fv);

nd = 2;
ncu = 3;    % num components of U
nch = ncu;    % num components of UHAT
nc = nch*(nd+1);    % num components of UDG

for n=1:1
    if n==1
        f = f1;
        filename1 = ['flux' num2str(nd) 'd' '.m'];
    elseif n==2
        f = f2;
        filename1 = ['flux' num2str(nd) 'd2' '.m'];
    end
    
    %%% compute Jacobian
    jac_f  = jacobian(f,udg);

    %%% And patch with vector zero or one to have the right sizes
    for ii = 1:size(f,1)
        for jj = 1:size(f,2)
            temp = ismember(symvar(f(ii,jj)),udg);
            if f(ii,jj)==0, f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, f(ii,jj) = (f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

    for ii = 1:size(jac_f,1)
        for jj = 1:size(jac_f,2)
            temp = ismember(symvar(jac_f(ii,jj)),udg);
            if jac_f(ii,jj)==0, jac_f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, jac_f(ii,jj) = (jac_f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end    
    
    %matlabFunction(f(1),'file','tmp.m','vars',{pg,udg,param,time,[zero one]},'outputs', {'f'});
    
    % generate a temporary matlab file      % NOTE: "p" was changed from the template
    matlabFunction(f(:),jac_f(:),'file','tmp.m','vars',{p,udg,param,time,[zero one]},'outputs', {'f','f_udg'});

    % open the file and modify it
    fid = fopen('tmp.m','r');
    gid = fopen(filename1,'wt');

    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        str = strrep(tline, 'tmp', strrep(filename1,'.m',''));
        str = strrep(str, 'TMP', upper(strrep(filename1,'.m','')));
        str = strrep(str, 'in1', 'p');          % NOTE: This was changed from the template
        str = strrep(str, 'IN1', 'P');        
        str = strrep(str, 'in2', 'udg');
        str = strrep(str, 'IN2', 'UDG');        
        str = strrep(str, 'in3', 'param');
        str = strrep(str, 'IN3', 'PARAM');                
        str = strrep(str, ',in5)', ')');                
        str = strrep(str, ',IN5)', ')');    
        str = strrep(str, 'param(:,1)', 'param{1}'); 
        str = strrep(str, 'param(:,2)', 'param{2}'); 
        str = strrep(str, 'param(:,3)', 'param{3}'); 
        str = strrep(str, 'param(:,4)', 'param{4}'); 
        str = strrep(str, 'param(:,5)', 'param{5}'); 
        str = strrep(str, 'param(:,6)', 'param{6}'); 
        str = strrep(str, 'param(:,7)', 'param{7}'); 
        str = strrep(str, 'param(:,8)', 'param{8}'); 
        str = strrep(str, 'param(:,9)', 'param{9}');     
        str = strrep(str, 'param(:,10)', 'param{10}');      % NOTE: These lines were added
        str = strrep(str, 'param(:,11)', 'param{11}');     
        str = strrep(str, 'param(:,12)', 'param{12}');     
        str = strrep(str, 'param(:,13)', 'param{13}');     
        str = strrep(str, 'param(:,14)', 'param{14}');     
        str = strrep(str, 'param(:,15)', 'param{15}');     
        str = strrep(str, 'param(:,16)', 'param{16}');     
        str = strrep(str, 'param(:,17)', 'param{17}');     
        str = strrep(str, 'param(:,18)', 'param{18}');     
        str = strrep(str, 'param(:,19)', 'param{19}');     
        str = strrep(str, 'in5(:,1)', 'zeros(ng,1)');
        str = strrep(str, 'in5(:,2)', 'ones(ng,1)');        
        if i==7
            str = '[ng,nc] = size(udg);';
            fprintf(gid, '%s\n', str);                  
            str = ['nch = ' num2str(nch) ';'];
            fprintf(gid, '%s\n', str);                  
            str = ['nd = ' num2str(nd) ';'];
        end
        fprintf(gid, '%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;
        %disp(str)
    end

    % NOTE: This will produce a function that has errors. You need to move the following "reshape" lies above the second "end" statement in order to be enclosed in the the function
    str = 'f = reshape(f,ng,nch,nd);';
    fprintf(gid, '%s\n', str);                  
    str = 'f_udg = reshape(f_udg,ng,nch,nd,nc);';
    fprintf(gid, '%s\n', str);                  

    fclose(fid);
    fclose(gid);
    delete('tmp.m');
end