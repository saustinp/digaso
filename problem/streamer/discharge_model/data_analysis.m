% time = 1ns
snapshot = 600;
cwi_slice = 3;

fname = sprintf("run_9-13-23/time%d.mat", snapshot);
UDG_snapshot = open(fname);
UDG = UDG_snapshot.UDG;
CWI_data = load('cwi_data/cwi_data_concatenated.mat');

porder = 2;
hybrid = 'hdg';
mesh = mkmesh_streamer_gmsh(porder, "streamer_16k-2.msh");
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% scaplot E field magnitude
% Er = UDG(:,6,:);
% Ez = UDG(:,9,:);
% normE = sqrt(Er.^2+Ez.^2)*3e6;
% figure(1); clf; scaplot(mesh,normE,[],2,0); axis equal; axis tight; axis on; colormap jet; title('|E|');
% 
% % scaplot ne
% figure(2); clf; scaplot(mesh,UDG(:,6,:)*1e12,[],2,0); axis equal; axis tight; axis on; colormap jet; title('ne');

% plot E field and ne

[bdry_pts, ne] = plot_bdry(mesh, UDG, 4, 1);
[~, Er_bdry] = plot_bdry(mesh, UDG, 4, 6);
[~, Ez_bdry] = plot_bdry(mesh, UDG, 4, 9);

normE_bdry = sqrt(Er_bdry.^2+Ez_bdry.^2)*3;

CWI_ycoord = CWI_data.cwi_data(:,2,cwi_slice);
CWI_bdry_E = CWI_data.cwi_data(:,6,cwi_slice);
CWI_ne = CWI_data.cwi_data(:,4,cwi_slice);

% figure(3);
% plot(CWI_ycoord*1e2, CWI_bdry_E/1e6, 'm', LineWidth=1.5, DisplayName='CWI, t=2ns'); hold on;
% plot(bdry_pts(:,2)*1e-2, normE_bdry, '--k', LineWidth=1.5, DisplayName='HDG, t=2ns'); hold on;
% legend();
% title('On-axis E field at t=0,1,2 ns');
% xlabel('z [cm]');
% ylabel('|E| [MV/m]');
% set(gca,'FontSize',16);


figure(4);
plot(CWI_ycoord*1e2, CWI_ne, 'm', LineWidth=1.5, DisplayName='CWI, t=2ns'); hold on;
plot(bdry_pts(:,2)*1e-2, ne*1e12, '--k', LineWidth=1.5, DisplayName='HDG, t=2ns'); hold on;
legend();
title('N_e at t=1,2 ns');
xlabel('z [cm]');
ylabel('n_e [m^{-3}]');
set(gca,'FontSize',16);