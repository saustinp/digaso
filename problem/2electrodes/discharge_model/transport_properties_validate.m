% Validate the transport properties

param = phys_param();     % Physics param loaded in a separate script
l_ref = param{1};
mu_ref = param{2};
E_ref = param{3};
e_eps0 = param{4};
normE = linspace(0,35,5000).*1e6;

% Compute quantities nondimensionally. Synthetic but it's important to check the functionality in-situ.
normE_tilde = normE./E_ref;
normE = normE_tilde.*E_ref;

mue_tilde = (2.3987.*normE.^(-.26))./mu_ref;
De_tilde = 4.3628e-3.*normE.^.22 ./ (l_ref.*mu_ref.*E_ref);
alpha = (1.1944e6+ 4.3666e26./normE.^3).*exp(-2.73e7./normE);
alpha_tile = alpha.*l_ref;

% Convert back to dimensional for plotting
mue = mue_tilde.*mu_ref;
De = De_tilde.*(l_ref.*mu_ref.*E_ref);
alpha = alpha_tile./l_ref;

% Plotting
close all;
figure();
plot(cwimue(:,1), cwimue(:,2), 'blue', LineWidth=2, DisplayName='CWI mu_e'); hold on;
plot(normE, mue, '--r', LineWidth=2, DisplayName='Curve fit');

legend();
title('Mobility coefficient');
xlabel('|E|, MV/m')
ylabel('mu_e, [m2/Vs]')

figure();
plot(cwidiff(:,1), cwidiff(:,2), 'blue', LineWidth=2, DisplayName='CWI D_e'); hold on;
plot(normE, De, '--r', LineWidth=2, DisplayName='Curve fit');

legend();
title('Diffusion coefficient');
xlabel('|E|, MV/m')
ylabel('D_e, [m2/s]')

figure();
plot(cwialpha(:,1), cwialpha(:,2), 'blue', LineWidth=2, DisplayName='CWI alpha'); hold on;
plot(normE, alpha, '--r', LineWidth=2, DisplayName='Curve fit');

legend();
title('Ionization coefficient');
xlabel('|E|, MV/m')
ylabel('alpha, [1/m]')