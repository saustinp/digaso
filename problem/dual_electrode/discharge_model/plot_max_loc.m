% Note: have to load the OLD mesh for this to work: 'streamer_mesh-16k.msh'

% mesh = mkmesh_streamer_gmsh(2, "streamer_16k.msh");
% old_max_loc_vec = zeros(76,1);
% old_max_vec = zeros(76,1);
% 
% for i=1:76
%     [~, Er] = plot_bdry(mesh, UDG_history(:,:,:,i), 4, 6);
%     [unique_pts, Ez] = plot_bdry(mesh, UDG_history(:,:,:,i), 4, 9);
% 
%     normE2d = sqrt(Er.^2+Ez.^2)*3e6;
%     [M,I] = max(normE2d);
%     old_max_vec(i) = M;
%     old_max_loc_vec(i) = unique_pts(I,2)*1e-4;
% 
% end
% 
% mesh = mkmesh_streamer_gmsh(2, "streamer_16k-2.msh");
% new_max_loc_vec = zeros(955,1);
% new_max_vec = zeros(955,1);
% for i=0:954
%     disp(i)
%     fname = sprintf('run_9-13-23/time%d.mat', i);
%     load(fname, "UDG");
%     [~, Er] = plot_bdry(mesh, UDG(:,:,:), 4, 6);
%     [unique_pts, Ez] = plot_bdry(mesh, UDG(:,:,:), 4, 9);
% 
%     normE2d = sqrt(Er.^2+Ez.^2)*3e6;
%     [M,I] = max(normE2d);
%     new_max_vec(i+1) = M;
%     new_max_loc_vec(i+1) = unique_pts(I,2)*1e-4;
%     clear("UDG")
% end

% plot_max_loc_new;       % Calls the script to process the most recent (9.26.23) run

figure(1); clf;
plot(cwi_summary(1:17,9)*100, cwi_summary(1:17,7), 'b', LineWidth=1.5, DisplayName='CWI/Truth'); hold on;
plot(old_max_loc_vec*100, old_max_vec, 'r', LineWidth=1.5, DisplayName='Mesh 1: Coarse 16k, 9/13/23');
plot(new_max_loc_vec*100, new_max_vec, 'm', LineWidth=1.5, DisplayName='Mesh 2: Coarse 16k, refined near IC, 9/14/23');
plot(max_loc_vec_2*100, max_vec_2, 'k', LineWidth=1.5, DisplayName='Mesh 2: Coarse 16k_adaptive, 9/26-27/23');

xlim([0,1])
xlabel('z [cm]');
title('On-axis E field max value and location');
ylabel('|E|, MV/m')
legend(Location='northwest')
set(gca,'FontSize',16);


