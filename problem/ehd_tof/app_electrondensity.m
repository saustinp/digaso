porder = 4;
hybrid = 'hdg';
nstage = 2; 
torder = 2;
tau = 1;
clear app;
app.axisymmetry = 1;

app.source = 'source2d';
app.flux = 'flux2d';
app.fbou = 'fbou_electrondensity';
app.fhat = 'fhat_electrondensity';
app.localsolve=1;

% Boundaries
% 1 Bottom electrode
% 2 Right farfield
% 3 Top electrode
% 4 Axisymmetry

% BC types
% 1 Axisymmetry
% 2 Homogeneous Neumann
% 3 Homogeneous Dirichlet

app.bcm = [2; 3; 2; 1;];
app.bcs = [0; 0; 0; 0];

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
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = app.nch;

ntime  = 68;
dt = 0.01*ones(ntime,1);

app.time = [];
app.dtfc = [];
app.alpha = [];

% Initializing data structures
app.arg = init_phys_param();     % Physics param loaded in a separate script
app.arg{end+1} = tau;

% mesh = mkmesh_tof(porder, "tof_mesh20k.msh");
mesh = mkmesh_rect2(11,41,porder,0,[0 1 0 2],0,1);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% Number density initialized to the same gaussian for both electrons and positives, and 0 for negatives.
% These have to be initialized in the same order that UDG is in
%                 ne_0,         q_ne_r0,        q_ne_z0
initu_func_set = {@initu_func; @initq_func_r; @initq_func_z;};
UDG = initu(mesh,initu_func_set,app.arg); 
UDG_history = zeros([size(UDG),ntime+1]);
UDG_history(:,:,:,1) = UDG;     % Adding in IC to the first snapshot
itime_restart = 0;

% Plot IC
% ne = UDG_history(:,1,:,1);
% scaplot(mesh,ne,[],0,1); axis equal; axis tight; colormap jet; title('u0');
% return;

UH=inituhat(master,mesh.elcon,UDG,app.ncu);

time = 0;
disp('Starting sim...')

for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime+itime_restart);

    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);
    time = time + dt(itime);

    UDG_history(:,:,:,itime+1+itime_restart) = UDG;
    figure(1); clf; scaplot(mesh,UDG(:,1,:),[],1,1); axis equal; axis tight; axis on;  
    
end
