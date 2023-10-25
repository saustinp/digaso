function param = init_phys_param()
    % Returns all of the simulation physical parameters that don't depend on the electric field

    De = 0.12; % m^2/s
    vd = 1.7e5; % m/s
    normE = 3.55e6; % V/m
    net_alpha = 5009.5;  % /m

    r_tip = 5e-4;            % [m] -- changed for this simulation!

    gamma = 0.001;           % Secondary electron emission coefficient [1/m]
    E_bd = 3e6;              % Breakdown E field in air [V/m]

    e = 1.6022e-19;          % Charge on electron [C]
    epsilon0 = 8.854e-12;    % absolute permittivity of air [C^2/(N*m^2)]
    n_ref = epsilon0*E_bd/(e*r_tip);  % Yes, this is redundant and could be recomputed from the above variables. But it saves having to recompute it each time in the functions.
    % ^=7.5357e+17
    t0 = 2e-9;   % s

    %        1    2    3         4       5      6     7      8    9
    param = {De, vd, normE, net_alpha, r_tip, gamma, E_bd, n_ref, t0};
end

% Reference list of physics parameters
% De = param(1);
% vd = param(2);
% normE = param(3);
% net_alpha = param(4);
% r_tip = param(5);
% gamma = param(6);
% E_bd = param(7);
% n_ref = param(8);
% t0 = param(9);

% Parameters to nondimensionalize
% Fundamental quantities: [L,T]
% ne: 1/L^3
% t: T
% vd: L/T
% De: L^2/T
% (alpha-eta): 1/L
% r,z: L

% Choose the domain radius for the reference legnth and the velocity/flow through time for the reference time

% Thus, the following nondimensional groups are formed: (_s = starred/nondimensional quantity)
% ne = ne_s/r^3
% t = t_s*(r/vd)
% De = De_s*r*vd
% (alpha-eta) = (alpha-eta)_s /r

% Thus, the PDE simplifies to:

% d ne_s/ dt_s   =  \nabla(ne_s - De_s \nabla(ne_s) ) = (alpha-eta)_s * ne_s
