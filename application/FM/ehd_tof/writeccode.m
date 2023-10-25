appname = 'ehd_tof';
nd=2;

syms x1 x2  % Variables stored in the p vector
syms u1 u2 u3     % Components of UDG
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10    % Components of the physics parameters. Last one is tau
syms nl1 nl2
syms uh1

param = [param1 param2 param3 param4 param5 param6 param7 param8 param9 param10];
udg = [u1 u2 u3];
pg = [x1 x2];

uh = uh1;
nl = [nl1 nl2];    

ncu = length(uh);
nc = length(udg);
nch = ncu;
                                       
% Read in values from the p vector
r = pg(1);
tau  = param(10).*r;

% Read in values from the u vector
ne = udg(1);
dne_dr = udg(2); % q is -grad(u)
dne_dz = udg(3);

% Loading physics parameters
De = param(1);
vd = param(2);
normE = param(3);
net_alpha = param(4);
r_tip = param(5);

De_s = De/(r_tip*vd);     % _s is "starred" or nondimensionalized quantity

fv = [De_s.*dne_dr, De_s.*dne_dz];
fi = [0, ne];

% Flux
f = r*(fi + fv);

% Fhat
fhi = [0, uh];
ff = r*(fhi + fv);
fhat = simplify(ff(1)*nl(1) + ff(2)*nl(2) + tau*(ne-uh));

% Source
s = r*net_alpha.*r_tip.*ne;



% filename1 = ['flux_' appname num2str(nd) 'd' '.c'];
% filename2 = ['source_' appname num2str(nd) 'd'  '.c'];
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