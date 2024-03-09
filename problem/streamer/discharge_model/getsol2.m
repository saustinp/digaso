% mesh = mkmesh_streamer_gmsh(2, "streamer_163k.msh");
% [UDG,UH] = getsolfrombinaryfile('./cluster_dl/soltime', './cluster_dl/streamersol',384,6,9,3,3,1,3600);   
% ne=UDG(:,1,:)*1e19;
% ni=UDG(:,2,:);
% Er = UDG(:,6,:);
% Ez = UDG(:,9,:);
% normE = sqrt(Er.^2+Ez.^2)*3e6;
% % 
% % % max(max(normE))
% % % max(max(ne))
% % figure(1); clf; colormap jet; scaplot(mesh,normE,[],0,0); axis equal; axis tight; title('ni-ne')
% % return;
% T = table2array(readtable('./cwi_data/comp_1e13_f_line_000013.txt'));
% % % 2:y, 3:ne, 6: |E|
% % 
% [unique_pts1, bdry_sol1] = plot_bdry(mesh, UDG, 4, 1);
% % [unique_pts, bdry_sol2] = plot_bdry(mesh, UDG, 4, 2);
% % [unique_pts, bdry_sol6] = plot_bdry(mesh, UDG, 4, 6);
% % [unique_pts, bdry_sol9] = plot_bdry(mesh, UDG, 4, 9);
% % % figure(2); 
% % 
% 
% mesh = mkmesh_streamer_gmsh(2, "streamer_350k.msh");
% [UDG,UH] = getsolfrombinaryfile('./cluster_dl350/soltime', './cluster_dl350/streamersol',384,6,9,3,3,1,3600);   
% ne=UDG(:,1,:)*1e19;
% ni=UDG(:,2,:);
% Er = UDG(:,6,:);
% Ez = UDG(:,9,:);
% normE = sqrt(Er.^2+Ez.^2)*3e6;
% [unique_pts2, bdry_sol2] = plot_bdry(mesh, UDG, 4, 1);



% load 'UDG3600_fullion.mat'
% UDG3600_fullion = UDG;
% [unique_pts, bdry_sol3] = plot_bdry(mesh, UDG3600_fullion, 4, 1);
% [unique_pts, bdry_sol4] = plot_bdry(mesh, UDG3600_fullion, 4, 2);
% % plot(unique_pts(:,2)*1e-2,sqrt(bdry_sol6.^2+bdry_sol9.^2), color='blue'); hold on;
% % plot(unique_pts(:,2)*1e-2,(bdry_sol2-bdry_sol1), color='blue'); hold on;
% plot(unique_pts(:,2)*1e-2,bdry_sol4-bdry_sol3, color='blue'); hold on;
% plot(unique_pts(:,2)*1e-2,bdry_sol2-bdry_sol1, color='red'); hold on;


% 
% 
plot(unique_pts1(:,2)*1e-2,bdry_sol1*1e19, color='blue'); hold on;
plot(unique_pts2(:,2)*1e-2,bdry_sol2*1e19, color='red', linewidth=2); hold on;
plot(T(:,2)*100, T(:,3), color='black', linewidth=2); hold on;
% plot(unique_pts(:,2)*1e-2, unique_pts(:,2)*0+1e13, color='black');
% xlim([0,.6])
% ylim([.5e13, 1.1e13])
ylabel('Ne')
xlabel('Distance, cm')
% 
legend('HDG :yuck:', 'HDG good', 'CWI/validate')
% 
% % Plot norm E, E source, to see how that damping effects the E field
% % irregularity. Also plot source term for the ion equation

% mesh = mkmesh_streamer_gmsh(porder, "streamer_163k.msh");
% load 'UDG3600_fullion.mat' UDG
% [unique_pts1, ~] = plot_bdry(mesh, UDG, 4, 1);
% diffy1=diff(unique_pts1(:,2));
% 

% return;
% mesh = mkmesh_streamer_gmsh(porder, "streamer_16k_fixed.msh");
% load './run030724_mat/time900.mat' UDG
% [unique_pts2, bdry_sol1] = plot_bdry(mesh, UDG, 4, 1);
% diffy2=diff(unique_pts2(:,2));
% % 
% % plot(unique_pts1(1:end-1,2), diffy1); hold on;
% % plot(unique_pts2(1:end-1,2), diffy2);
% % 
% figure();
% plot(unique_pts2(:,2)*1e-2,bdry_sol1, color='red'); hold on;
% scatter(unique_pts2(:,2)*1e-2,bdry_sol1, color='blue'); hold on;
% plot(unique_pts2(1:end-1,2)*1e-2,diffy2*500, color='black'); hold on;