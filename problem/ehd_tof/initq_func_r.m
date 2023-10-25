function q0 = initq_func_r(dgnodes, param)
    r = dgnodes(:,1,:);
    z = dgnodes(:,2,:);

    % Physics parameters
    De = param{1};
    vd = param{2};
    normE = param{3};
    net_alpha = param{4};
    r_tip = param{5};
    gamma = param{6};
    E_bd = param{7};
    n_ref = param{8};
    t0 = param{9};

    t = t0./(r_tip./vd);

    q0 = (r.*r_tip.*vd.*exp(net_alpha.*r_tip.*t - (r_tip.*vd.*((t - z).^2 + r.^2))./(4.*De.*t)))./(2.*De.*t.*((4.*De.*pi.*t)./(r_tip.*vd)).^(3./2));
end