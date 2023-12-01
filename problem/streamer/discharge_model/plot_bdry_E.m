itime = 260;

UDG_snapshot = UDG_history(:,:,:,itime);



[~, ne] = plot_bdry(mesh, snapshot_230, 4, 1);
[~, ni] = plot_bdry(mesh, snapshot_230, 4, 2);

[~, Er] = plot_bdry(mesh, UDG_snapshot, 4, 6);
[unique_pts, Ez] = plot_bdry(mesh, UDG_snapshot, 4, 9);
normE2d = sqrt(Er.^2+Ez.^2);

figure(1);
% plot(unique_pts(:,2)*1e-2, normE2d*3e6, LineWidth=1.5, DisplayName='HDG'); hold on;
plot(unique_pts(:,2)*1e-2, normE2d*3e6, LineWidth=1.5); hold on;
% plot(unique_pts(:,2)*1e-2, ni-ne, LineWidth=1.5); hold on;

legend();
% title('On-axis E field');
xlabel('z [cm]')
% ylabel('|E|, MV/m')