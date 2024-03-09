function u0 = initu_func_electrons(dgnodes, param)
    r_tilde = dgnodes(:,1,:);
    z_tilde = dgnodes(:,2,:);

    % Physics parameters
    l_ref = param{1};
    N0 = param{6};
    sigma0 = param{8};
    n_background = param{9};

    z0e = -53.585*l_ref; % m
    N0_tilde = N0*(l_ref^3);
    n_background_tilde = n_background*(l_ref^3);
    z0e_tilde = z0e/l_ref;
    sigma0_tilde = sigma0/l_ref;

    u0 = n_background_tilde + N0_tilde * exp(-((z_tilde - z0e_tilde).^2 + r_tilde.^2)/(sigma0_tilde^2));
    % u0 = n_background_tilde;
end