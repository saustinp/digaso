function param = phys_param()
    % Returns all of the simulation physical parameters that don't depend on the electric field

    kappa = 1;  % Diffusivity
    l_ref = 1e-4;   % m
    mu_ref = 0.05;  % m^2/(V-s)
    E_ref = 3e6;    % V/m
    e_eps0 = 1.80955e-8;    % Quantity e/epsilon0
    phi0 = 5e3;     % V
    N0 = 3e18; % 1/m^3
    z0 = -3.58*l_ref; % m
    sigma0 = 4e-4;  % m
    z0e = -53.585*l_ref; % m

    %          1     2      3       4      5       6    7    8    9      10
    param = {kappa, l_ref, mu_ref, E_ref, e_eps0, phi0, N0, z0, sigma0, z0e};
end

% Reference list of physics parameters
% kappa = param{1};
% l_ref = param{2};
% mu_ref = param{3};
% E_ref = param{4};
% e_eps0 = param{5};
% phi0 = param{6};
% N0 = param{7};
% z0 = param{8};
% sigma0 = param{9};