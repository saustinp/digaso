syms n_background_tilde N0_tilde z_tilde z0_tilde r_tilde sigma0_tilde

% First, write the nondimensional quantities here:

u0 = n_background_tilde + N0_tilde * exp(-((z_tilde - z0_tilde).^2 + r_tilde.^2)/(sigma0_tilde^2));

du0_dr = -diff(u0, r_tilde)   % Need the negative for q. The result is then printed out to the terminal and copy/pasted into the setup script
du0_dz = -diff(u0, z_tilde)

% Then:
% 1. copy/paste the output below
% 2. Copy/paste together into an empty file and find/replace *, /, ^ with elementwise counterparts
% 3. Paste into the respective q0 scripts

% du0_dr = (2.*N0_tilde.*r_tilde.*exp(-((z0_tilde - z_tilde).^2 + r_tilde.^2)./sigma0_tilde.^2))./sigma0_tilde.^2
% du0_dz = -(N0_tilde.*exp(-((z0_tilde - z_tilde).^2 + r_tilde.^2)./sigma0_tilde.^2).*(2.*z0_tilde - 2.*z_tilde))./sigma0_tilde.^2