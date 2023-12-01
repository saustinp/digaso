setapplicationpath('FM/ehd_tof');

fileName = 'tof';

porder = 4;
nstage = 2; 
torder = 2;
hybrid = 'hdg';
app.porder=porder;
app.torder=torder;
app.nstage=nstage;

% Physics params
% De = 0.12; % m^2/s
% vd = 1.7e5; % m/s
% normE = 3.55e6; % V/m
% net_alpha = 5009.5;  % /m
% r_tip = 5e-4;            % [m] -- changed for this simulation!
% gamma = 0.001;           % Secondary electron emission coefficient [1/m]
% E_bd = 3e6;              % Breakdown E field in air [V/m]
% e = 1.6022e-19;          % Charge on electron [C]
% epsilon0 = 8.854e-12;    % absolute permittivity of air [C^2/(N*m^2)]
% n_ref = epsilon0*E_bd/(e*r_tip);  % Yes, this is redundant and could be recomputed from the above variables. But it saves having to recompute it each time in the functions.
% % ^=7.5357e+17
% t0 = 2e-9;   % s
% tau = 1;
% %        1    2    3         4       5      6     7      8    9    10
% app.arg = {De, vd, normE, net_alpha, r_tip, gamma, E_bd, n_ref, t0, tau};

app.ncu = 3;
app.nd = 2;
app.arg = init_phys_param();     % Physics param loaded in a separate script
app.arg{end+1} = 1;

app.appname = 'ehd_tof';
app.axisymmetry = 1;

app.localsolve=1;

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


app.iterative = 0;      % This flag is in both the rans/naca0012 and the poisson case
app.uqpk    = 0;
app.morder          = [porder porder];
app.porder          = [porder porder];
nodetype    = 1;
app.nodetype        = nodetype;
app.pgauss          = 2*[porder porder];
app.pgaussR         = 2*[porder porder];
app.overlappinglevel= 1;
app.preconditioner  = 0;
app.ncq     = app.nd*app.ncu; 
app.ncp     = 0; 
app.nco     = 0;    % Additional field for defining AV
app.quadtype        = [0 0];


app.nd = 2;
app.nch  = 1;                       % Number of componets of UH
app.ncu = app.nch;
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG

app.time = 0;
ntime  = 68;    % Number of timesteps to take
app.dt = 0.01*ones(ntime,1);

% Initializing data structures
% app.arg = init_phys_param();     % Physics param loaded in a separate script
% app.arg{end+1} = tau;

% mesh = mkmesh_tof(porder, "tof_mesh20k.msh");
mesh = mkmesh_rect2(11,41,porder,0,[0 1 0 2],0,1);
app.ncd     = size(mesh.dgnodes,2);
elementtype = mesh.elemtype*ones(mesh.ne,1);
check = 0;

master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% Number density initialized to the same gaussian for both electrons and positives, and 0 for negatives.
% These have to be initialized in the same order that UDG is in
%                 ne_0,         q_ne_r0,        q_ne_z0
initu_func_set = {@initu_func; @initq_func_r; @initq_func_z;};
UDG0 = initu(mesh,initu_func_set,app.arg); 
UH0=inituhat(master,mesh.elcon,UDG0,app.ncu);

% [Q,~,~] = getq(master,mesh,UDG0,UH0);
% UDG0 = cat(2,UDG0,Q);

% ------- GMRES Parameters --------- %
app.restart   = 200;
app.gmrestol  = 1e-12;
app.gmresiter = 2000;

% ------- Newton Parameters ---------- %
app.newtoniter = 100;  % def 10
app.newtontol  = 1e-8; % def 1e-7

% ------- Number of processors ------- %
nproc       = 1;
app.nfile   = nproc;

delete *.bin
if nproc>1 % parallel preprocessing
    apppar = digasopre(app,'ehd_tof',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,mesh.bndexpr,[],nproc,0,check,mesh.perm);
    apppar.fileout = 'ehd_tof';
    return;
else  % serial preprocessing
    appser = digasopre(app,'ehd_tof',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,mesh.bndexpr,[],nproc,0,check,mesh.perm);
    appser.fileout = 'ehd_tof';
    return;
end

return;


