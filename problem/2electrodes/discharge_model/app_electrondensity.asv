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

app.bcm = [1; 2; 1; 3;];            % Mod for testing the C++ code
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

% app.check_volint_flg = 1;
% app.check_fhat_flg = 1;
% app.check_fbou_flg = 4;
app.debug_digaso = 1;

ntime  = 2000;
dt = 5e-3*ones(ntime,1);      % With choice of nondimensionalization, t_tilde=1 => 6.67e-10s

app.time = [];
app.dtfc = [];
app.alpha = [];

% mesh = mkmesh_rect(41,81,porder,0,[0 125 0 125],0,1);
mesh = mkmesh_streamer_gmsh(porder, "streamer_16k_fixed.msh");
linearmesh = mkmesh_streamer_gmsh(1, "streamer_16k_fixed.msh");

% mesh = mkmesh_streamer_gmsh(2, "streamer_16k-3.msh");
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid,app);
% meshplot(mesh);
% return;

% Uncomment to start sim from 0
% Number density initialized to the same gaussian for both electrons and positives, and 0 for negatives.
% These have to be initialized in the same order that UDG is in
%                  ne_0,                     np_0,   phi_0,  q_ne_r0,q_np_r0,Er0,       q_ne_z0,q_np_z0,Ez0
% 
initu_func_set = {@initu_func_electrons;@initu_func_ions;0;   0;@initq_func_ions_r;0;    0;@initq_func_ions_z;0};

load '../poissonIC.mat';
UDG_poisson = UDG;      % Load in the poisson as UDG
UDG = initu(mesh,initu_func_set,app.arg);
UDG(:,[3,6,9],:) = UDG_poisson;
UH=inituhat(master,mesh.elcon,UDG,app.ncu);
itime_restart = 0;

% Restart
% UDG=load("restartSol.mat");
% load restartSol.mat
% itime_restart = 501;
% UH=inituhat(master,mesh.elcon,UDG,app.ncu);

% diary run_11_7_23/run_11_7_23.txt;
% diary on;
time = itime_restart*dt(1);     % Need to change if non-constant dt
disp('Starting sim...')

% save 'run_11_7_23/time501' UDG;
% return;
for itime = 1:ntime
    diary on;
    fprintf('Timestep :  %d\n', itime+itime_restart);

    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime+itime_restart),nstage,torder,linearmesh);
    time = time + dt(itime+itime_restart);
    % fname_out = 'run_11_7_23/time' + string(itime+itime_restart);
    % save(fname_out, "UDG");

    disp('ne max')
    disp(max(max(UDG(:,1,:))))
    disp('E max');

    Er = UDG(:,6,:);
    Ez = UDG(:,9,:);
    normE = sqrt(Er.^2+Ez.^2)*3e6;
    % scaplot(mesh,normE,[],0,0); axis equal; axis tight; colormap jet; title('|E|');
    disp(max(max(normE)))
    % diary off;

end
