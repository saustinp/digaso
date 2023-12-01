appname = 'streamer';
nd=2;

syms x1 x2  % Variables stored in the p vector
syms u1 u2 u3 u4 u5 u6 u7 u8 u9     % Components of UDG
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10    % Components of the physics parameters
syms nl1 nl2
syms uh1 uh2 uh3

pg = [x1 x2];
udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9];
param = [param1 param2 param3 param4 param5 param6 param7 param8 param9 param10];
nl = [nl1 nl2];    
uh = [uh1 uh2 uh3];

ncu = length(uh);
nc = length(udg);
nch = ncu;
                                       
% Read in values from the p vector
r_tilde = pg(1);
tau  = param(10);

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
e_eps0 = param(4);

% Compute transport and source coefficients
normE_tilde = sqrt(Er_tilde^2 + Ez_tilde^2);
normE = normE_tilde*E_ref;

mue_tilde = (2.3987*normE^(-.26))/mu_ref;
De_tilde = 4.3628e-3*normE^.22 / (l_ref*mu_ref*E_ref);
alpha = (1.1944e6+ 4.3666e26/normE^3)*exp(-2.73e7/normE);
alpha_tile = alpha*l_ref;
eta_tilde = 340.75*l_ref;
alpha_bar_tilde = alpha_tile-eta_tilde;

% Flux
fv = [De_tilde.*dne_dr_tilde,0,Er_tilde,...
      De_tilde.*dne_dz_tilde,0,Ez_tilde];       % The negative sign is included in eqns 1 and 3 bc q=-grad(u)

fi = [-Er_tilde*mue_tilde*ne_tilde,0,0,...
      -Ez_tilde*mue_tilde*ne_tilde,0,0];

f = r_tilde*(fi + fv);

% Fhat
fhi = [-Er_tilde*mue_tilde*uh(1),0,0,...
       -Ez_tilde*mue_tilde*uh(1),0,0];
       
ff = r_tilde*(fhi + fv);

fhat = [simplify(ff(1)*nl(1) + ff(4)*nl(2) + tau*(ne_tilde-uh(1))),...
        simplify(ff(2)*nl(1) + ff(5)*nl(2) + tau*(ni_tilde-uh(2))),...
        simplify(ff(3)*nl(1) + ff(6)*nl(2) + tau*(phi_tilde-uh(3)))];

% Source
charge_prefactor_tilde = e_eps0/(E_ref*l_ref^2);
se = alpha_bar_tilde*mue_tilde*normE_tilde*ne_tilde;
sphi = charge_prefactor_tilde*(ni_tilde - ne_tilde);    % Note sign flip in source term because the negative sign is absorbed in the diffusive flux in the LHS

s = r_tilde*[se se sphi];

filename1 = ['fluxG_' appname num2str(nd) 'd' '.c'];
filename2 = ['sourceG_' appname num2str(nd) 'd'  '.c'];
filename3 = ['fhat_' appname num2str(nd) 'd' '.c'];

genccode; % generate source codes

% filename1 = ['fluxonly_' appname num2str(nd) 'd' '.c'];
% filename2 = ['sourceonly_' appname num2str(nd) 'd'  '.c'];
% filename3 = ['fhatonly_' appname num2str(nd) 'd' '.c'];
% genccode_withoutjac; % generate source codes

% generate an application file
% gid = fopen('fluxes.c','w');
% fid = fopen('../../fluxes_template.c','r');    
% tline = fgetl(fid); 
% while ischar(tline)        
%     str = tline;
%     str = strrep(str, 'appname', appname);    
%     fprintf(gid, '%s\n', str);    
%     tline = fgetl(fid);            
% end
% fclose(fid);
% fclose(gid);