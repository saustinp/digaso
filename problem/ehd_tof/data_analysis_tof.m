% Data analysis script for the electron time of flight simulation

[unique_pts, bdry_sol1] = plot_bdry(mesh, UDG_history(:,:,:,1), 4, 1);

y_coord = unique_pts(:,2);
bdry_sol_all_time = zeros(size(bdry_sol1,1), ntime+1);
bdry_sol_all_time(:,1) = bdry_sol1;

for j=1:ntime
    [unique_pts, bdry_sol] = plot_bdry(mesh, UDG_history(:,:,:,j+1), 4, 1);
    bdry_sol_all_time(:,j+1) = bdry_sol;
end

% Rescale to dimensional units
lref = 5e-4;
De = 0.12; % m^2/s
vd = 1.7e5; % m/s
net_alpha = 5009.5;  % /m

t_ref = lref/vd;
t0 = 2e-9;  % in s

dt_s = dt*t_ref;
bdry_sol_all_time = bdry_sol_all_time*(lref^-3);
t_vec = t0 + linspace(0, dt(1)*ntime, ntime+1)*t_ref;     % Adding in one extra point (and starting from 0) to account for the IC being added into the history vector
y_coord = y_coord*lref;

% Analytical solution
exact_sol_all_time = zeros(size(bdry_sol_all_time));
for i=1:size(t_vec,2)
    t = t_vec(i);
    exact_sol_all_time(:,i) = (4*pi*De*t)^(-3/2) * exp(-(y_coord-vd*t).^2/(4*De*t) + net_alpha*vd*t);
end

% h(1) = plot(y_coord, exact_sol_all_time(:,1), 'b', 'DisplayName', 'Exact Solution'); hold on;
% pause;
% h(2) = plot(y_coord, bdry_sol_all_time(:,1), '--r', 'DisplayName', 'HDG'); hold on;
% pause;
% 
% plot(y_coord, exact_sol_all_time(:,2), 'b'); hold on;
% pause;
% plot(y_coord, bdry_sol_all_time(:,2), '--r'); hold on;
% pause;
% 
% plot(y_coord, exact_sol_all_time(:,3), 'b'); hold on;
% pause;
% plot(y_coord, bdry_sol_all_time(:,3), '--r'); hold on;
% pause;
% 
% legend([h(1), h(2)], 'Exact Solution', 'HDG')
% return;

y_coord = y_coord*1e3;      % Conversion to mm for plotting
figure();

lw=1.5;
h(1) = plot(y_coord, exact_sol_all_time(:,1), 'black', 'DisplayName', 'Exact Solution', LineWidth=lw); hold on;
% h(2) = plot(y_coord, bdry_sol_all_time(:,1), '--r', 'DisplayName', 'HDG', LineWidth=lw); hold on;

h(3) = plot(y_coord, exact_sol_all_time(:,18), 'b', 'DisplayName', 'Exact Solution', LineWidth=lw); hold on;
h(4) = plot(y_coord, bdry_sol_all_time(:,18), '--r', 'DisplayName', 'HDG', LineWidth=lw); hold on;

h(5) = plot(y_coord, exact_sol_all_time(:,35), 'b', 'DisplayName', 'Exact Solution', LineWidth=lw); hold on;
h(6) = plot(y_coord, bdry_sol_all_time(:,35), '--r', 'DisplayName', 'HDG', LineWidth=lw); hold on;

h(7) = plot(y_coord, exact_sol_all_time(:,52), 'b', 'DisplayName', 'Exact Solution', LineWidth=lw); hold on;
h(8) = plot(y_coord, bdry_sol_all_time(:,52), '--r', 'DisplayName', 'HDG', LineWidth=lw); hold on;

h(9) = plot(y_coord, exact_sol_all_time(:,69), 'b', 'DisplayName', 'Exact Solution', LineWidth=lw); hold on;
h(10) = plot(y_coord, bdry_sol_all_time(:,69), '--r', 'DisplayName', 'HDG', LineWidth=lw); hold on;

ylabel('{\it n_e} [10^{13} m^{-3}]');
xlabel('z [mm]');
set(gca,'fontsize', 14);

ylim([0, 7e13]);

max1 = max(exact_sol_all_time(:,1));
max17 = max(exact_sol_all_time(:,18));
max34 = max(exact_sol_all_time(:,35));
max51 = max(exact_sol_all_time(:,52));
max68 = max(exact_sol_all_time(:,69));

xoffset = 0.025;
yoffset = .25e13;
text(.34-xoffset, max1+yoffset, '2 ns')
text(.425-xoffset, max17+yoffset, '2.5 ns')
text(.51-xoffset, max34+yoffset, '3 ns')
text(.595-xoffset, max51+yoffset, '3.5 ns')
text(.68-xoffset, max68+yoffset, '4 ns')

legend([h(3), h(4)], 'Exact Solution', 'HDG', Location='northwest');

