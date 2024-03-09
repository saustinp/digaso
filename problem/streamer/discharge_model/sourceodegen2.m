% These inputs must match the dimensions of the input vectors in the final output function. For example:
% function [f,f_udg] = flux_electrondensity(p,udg,param,time)

syms x1 x2  % Variables stored in the p vector
syms u1 u2 u3 u4 u5 u6     % Components of UDG
syms o1 o2                   % Components of ODG
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10 param11   % Components of the physics parameters
syms factor1 factor2
syms zero one

param = [param1 param2 param3 param4 param5 param6 param7 param8 param9 param10 param11];
factor = [factor1 factor2];
udg = [u1 u2 u3 u4 u5 u6];
odg = [o1 o2];
p = [x1 x2];

% Read in values from the p vector
r_tilde = p(1);

% Read in values from the u vector
ne_tilde      = udg(1);
phi_tilde     = udg(2);
dne_dr_tilde  = udg(3);
Er_tilde      = udg(4);
dne_dz_tilde  = udg(5);
Ez_tilde      = udg(6);

% Load physics param
l_ref = param(1);
mu_ref = param(2);
E_ref = param(3);
e_eps0 = param(4);
ne_star = param(10);

% Compute transport coefficients
normE_tilde = sqrt(Er_tilde^2 + Ez_tilde^2);
normE = normE_tilde*E_ref;

mue_tilde = (2.3987*normE^(-.26))/mu_ref;
De_tilde = 4.3628e-3*normE^.22 / (l_ref*mu_ref*E_ref);
alpha = (1.1944e6+ 4.3666e26/normE^3)*exp(-2.73e7/normE);
alpha_tile = alpha*l_ref;
eta_tilde = 340.75*l_ref;
alpha_bar_tilde = alpha_tile-eta_tilde;

charge_prefactor_tilde = e_eps0/(E_ref*l_ref^2);

se = alpha_bar_tilde*mue_tilde*normE_tilde*ne_tilde;
f1 = se;

nd = 2;
ncu = 2;    % num components of U
nch = 1;    % num components of UHAT
nc = ncu*(nd+1);    % num components of UDG

for n=1:1
    if n==1
        f = f1;
        filename1 = ['sourceode' num2str(nd) 'd2' '.m'];
    elseif n==2
        f = f2;
        filename1 = ['sourceode' num2str(nd) 'd2' '.m'];
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
    
    % generate a temporary matlab file         $ NOTE: changed "pg" -> "p"
    matlabFunction(f(:),jac_f(:),'file','tmp.m','vars',{p,udg,param,time,[zero one]},'outputs', {'s','s_udg'});          % NOTE: This was changed from the template

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
        str = strrep(str, 'param(:,10)', 'param{10}');     
        str = strrep(str, 'param(:,11)', 'param{11}');       % NOTE: These lines were added
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
    str = 's = reshape(s,ng,nch);';          % NOTE: This was changed from the template
    fprintf(gid, '%s\n', str);         
    
    fprintf(gid, '%s\n', 'if nargout > 1');                  
    str = 's_udg = reshape(s_udg,ng,nch,nc);';          % NOTE: This was changed from the template
    fprintf(gid, '%s\n', str);    
    fprintf(gid, '%s\n', 'end');                  

    fclose(fid);
    fclose(gid);
    delete('tmp.m');
end