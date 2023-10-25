function u0 = initu_func(dgnodes, param)
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
    mue_ref = param{10};

    t = t0/(r_tip/vd);

    u0 = (4.*pi.*De.*t./r_tip.^2.* (r_tip./vd)).^(-3./2) .* exp(((z-t).^2 + r.^2)./(-4.*De.*t./r_tip.^2.* r_tip./vd) + net_alpha.*vd.*t .* r_tip./vd);
end