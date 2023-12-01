function q0 = initq_func_ions_z(dgnodes, param)
    r_tilde = dgnodes(:,1,:);
    z_tilde = dgnodes(:,2,:);

    % Physics parameters
    l_ref = param{1};
    N0 = param{6};
    z0 = param{7};
    sigma0 = param{8};

    N0_tilde = N0*(l_ref^3);
    z0_tilde = z0/l_ref;
    sigma0_tilde = sigma0/l_ref;

    q0 = -(N0_tilde.*exp(-((z0_tilde - z_tilde).^2 + r_tilde.^2)./sigma0_tilde.^2).*(2.*z0_tilde - 2.*z_tilde))./sigma0_tilde.^2;
end