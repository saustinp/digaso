porder = 2;
hybrid = 'hdg';
nstage = 1; 
torder = 1;
tau = 1;
clear app;
app.axisymmetry = 1;

app.arg = phys_param();     % Physics param loaded in a separate script
app.arg{end+1} = tau;

app.source = 'source2d';
app.flux = 'flux2d';
app.fbou = 'fbou_electrondensity';
app.fhat = 'fhat_axisymmetric';
app.localsolve=1;

% Boundaries
% 1 Bottom electrode
% 2 Right farfield
% 3 Top electrode
% 4 Symmetry

% BC types
% 1 Symmetry
% 2 Homogeneous Neumann
% 3 Homogeneous Dirichlet

l_ref = app.arg{1};
E_ref = app.arg{3};
phi0 = app.arg{5};
phi0_tilde = phi0/(E_ref*l_ref);

app.bcm = [1; 2; 1; 3;];
app.bcs = [0;0;phi0_tilde;0];
app.fcu_vector = [1;1;0];
% app.bcs = [0;0;0;0];

app.bcd = [];
app.bcv = [];

app.hybrid = hybrid;
app.tdep = true;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.fc_q = 1;          
app.fc_u = 1;
app.fc_p = 0;

app.nd = 2;
app.nch  = 3;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = app.nch;

ntime  = 2000;
dt = 5e-3*ones(ntime,1);      % With choice of nondimensionalization, t_tilde=1 => 6.67e-10s

app.time = [];
app.dtfc = [];
app.alpha = [];

% mesh = mkmesh_rect(41,81,porder,0,[0 125 0 125],0,1);
mesh = mkmesh_streamer_gmsh(porder, "streamer_16k-2.msh");
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
% meshplot(mesh);
% return;

% Uncomment to start sim from 0
% Number density initialized to the same gaussian for both electrons and positives, and 0 for negatives.
% These have to be initialized in the same order that UDG is in
%                  ne_0,                     np_0,   phi_0,  q_ne_r0,q_np_r0,Er0,       q_ne_z0,q_np_z0,Ez0

initu_func_set = {@initu_func_electrons;@initu_func_ions;0;   0;@initq_func_ions_r;0;    0;@initq_func_ions_z;0};

load '../poissonIC.mat';
UDG_poisson = UDG;      % Load in the poisson as UDG
UDG = initu(mesh,initu_func_set,app.arg);
UDG(:,[3,6,9],:) = UDG_poisson;
UH=inituhat(master,mesh.elcon,UDG,app.ncu);

% Initialize UDG_history array
% UDG_history = zeros([size(UDG),ntime+1]);
% UDG_history(:,:,:,1) = UDG;     % Adding in IC to the first snapshot
itime_restart = 0;

% Restart
% itime_restart = 35;
% UDG = UDG_history(:,:,:,itime_restart);

% Plot IC
% close all;
% figure();
% scaplot(mesh,UDG(:,1,:),[],0,1); axis equal; axis tight; colormap jet; title('n_{e,0}');
% figure();
% scaplot(mesh,UDG(:,2,:),[],0,0); axis equal; axis tight; colormap jet; title('n_{i,0}');
% figure();
% scaplot(mesh,UDG(:,3,:),[],0,0); axis equal; axis tight; colormap jet; title('phi0');
% figure();
% scaplot(mesh,UDG(:,4,:),[],0,0); axis equal; axis tight; colormap jet; title('dn_e/dr_0');
% figure();
% scaplot(mesh,UDG(:,5,:),[],0,0); axis equal; axis tight; colormap jet; title('dn_i/dr_0');
% figure();
% scaplot(mesh,UDG(:,6,:),[],0,0); axis equal; axis tight; colormap jet; title('Er_0');
% figure();
% scaplot(mesh,UDG(:,7,:),[],0,0); axis equal; axis tight; colormap jet; title('dn_e/dz_0');
% figure();
% scaplot(mesh,UDG(:,8,:),[],0,0); axis equal; axis tight; colormap jet; title('dn_i/dz_0');
% figure();
% scaplot(mesh,UDG(:,9,:),[],0,0); axis equal; axis tight; colormap jet; title('Ez_0');
% return;

diary run_10-13-23.txt;
diary on;
time = itime_restart*dt(1);     % Need to change if non-constant dt
disp('Starting sim...')

% save 'run_9-26-23/time0' UDG;

for itime = 1:ntime
    diary on;
    fprintf('Timestep :  %d\n', itime+itime_restart);

    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime+itime_restart),nstage,torder);
    time = time + dt(itime+itime_restart);
    fname_out = 'run_10-13-23/time' + string(itime+itime_restart);
    save(fname_out, "UDG");

    % UDG_history(:,:,:,itime+1+itime_restart) = UDG;
    % 
    % Er = UDG(:,6,:);
    % Ez = UDG(:,9,:);
    % normE = sqrt(Er.^2+Ez.^2)*3e6;
    % figure(1); clf; scaplot(mesh,normE,[],2,0); axis equal; axis tight; axis on; colormap jet;

    % figure(1); clf;    [~, bdry_sol11] = plot_bdry(mesh, UDG, 4, 6);
    % [unique_pts, bdry_sol12] = plot_bdry(mesh, UDG, 4, 9);
    % normE2d1 = sqrt(bdry_sol11.^2 + bdry_sol12.^2);
    % plot(unique_pts(:,2), normE2d1*3e6, DisplayName='current timestep'); hold on;
    % [~, bdry_sol01] = plot_bdry(mesh, UDG_poisson, 4, 2);
    % [~, bdry_sol02] = plot_bdry(mesh, UDG_poisson, 4, 3);
    % normE2d0 = sqrt(bdry_sol01.^2 + bdry_sol02.^2);
    % plot(unique_pts(:,2), normE2d0*3e6, DisplayName='T0'); hold on;
    % legend();

    disp('ne max')
    disp(max(max(UDG(:,1,:))))
    disp('E max');

    Er = UDG(:,6,:);
    Ez = UDG(:,9,:);
    normE = sqrt(Er.^2+Ez.^2)*3e6;
    scaplot(mesh,normE,[],0,0); axis equal; axis tight; colormap jet; title('|E|');
    disp(max(max(normE)))
    diary off;

end
