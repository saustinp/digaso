function [sr,sr_udg] = source_poisson(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

r_tilde = p(:,1);
z_tilde = p(:,2);

% Physics parameters
l_ref = param{2};
E_ref = param{4};
e_eps0 = param{5};
N0 = param{7};
z0 = param{8};
sigma0 = param{9};

N0_tilde = N0*(l_ref^3);
z0_tilde = z0/l_ref;
sigma0_tilde = sigma0/l_ref;

sr = r_tilde.*(e_eps0/(E_ref*l_ref^2)).*N0_tilde.*exp(-((z_tilde - z0_tilde).^2 + r_tilde.^2)/(sigma0_tilde^2));

% ni_discharge = p(:,3);
% ne_discharge = p(:,4);
% sr = r_tilde.*(e_eps0/(E_ref*l_ref^2)).*   (ni_discharge - ne_discharge);

sr_udg = zeros(ng,nch,nc);