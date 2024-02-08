[~, bdry_sol2] = plot_bdry(mesh, UDG, 4, 2);
[unique_pts, bdry_sol3] = plot_bdry(mesh, UDG, 4, 3);
normE2d = sqrt(bdry_sol2.^2 + bdry_sol3.^2);

[~, bdry_sol230] = plot_bdry(mesh, UDG_history(:,:,:,230), 4, 9);

figure();
plot(cwidata(:,2)*1e2, cwidata(:,6), 'blue', LineWidth=2, DisplayName='CWI/Truth'); xlim([0,1.25]); hold on;
% plot(unique_pts(:,2)*1e-2, bdry_sol3*3e6, LineWidth=2); hold on;
% plot(unique_pts(:,2)*1e-2, bdry_sol230*3e6, LineWidth=2);
plot(unique_pts(:,2)*1e-2, normE2d*3e6, '--r', LineWidth=2, DisplayName='HDG');

legend();
title('On-axis E field, t=0ns');
xlabel('z [cm]')
ylabel('|E|, MV/m')