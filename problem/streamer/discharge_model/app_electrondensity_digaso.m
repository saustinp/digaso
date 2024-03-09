setapplicationpath('FM/streamer');
fileName = 'streamer';

clear app;
porder = 2;
app.nstage = 2; 
app.torder = 2;
app.hybrid = 'hdg';

ne_star = 10^7; % this constant is used to normalize ne and ni
% ne_star = 1; % this constant is used to normalize ne and ni

% Init phys params
tau = 1;
% tau = 1e6;
app.arg = phys_param();     % Physics param loaded in a separate script
app.arg{end+1} = ne_star;
app.arg{end+1} = tau;

app.appname = 'streamer';
app.axisymmetry = 1;

app.localsolve=1;

l_ref = app.arg{1};
E_ref = app.arg{3};
phi0 = app.arg{5};
phi0_tilde = phi0/(E_ref*l_ref);

app.bcm = [1; 2; 1; 3;];
%app.bcs = [0;0;phi0_tilde;0];
app.bcs  = [[10/ne_star 10/ne_star 0]; [0 0 0]; [10/ne_star 10/ne_star phi0_tilde]; [0 0 0]];
app.fcu_vector = [1;1;0];
% app.bcs = [0;0;0;0];

app.bcd = [];
app.bcv = [];

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

app.time = 0;
ntime  = 6000;
app.dt = 5e-3*ones(ntime,1);      % With choice of nondimensionalization, t_tilde=1 => 6.67e-10s

% mesh = mkmesh_streamer_gmsh(porder, "streamer_16k_fixed.msh");
% mesh = mkmesh_streamer_gmsh(porder, "streamer_57k.msh");
% mesh = mkmesh_streamer_gmsh(porder, "streamer_163k.msh");
mesh = mkmesh_streamer_gmsh(porder, "streamer_350k.msh");

master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,app.hybrid);

% Number density initialized to the same gaussian for both electrons and positives, and 0 for negatives.
% These have to be initialized in the same order that UDG is in
%                  ne_0,                     np_0,     phi_0, q_ne_r0,q_np_r0,    Er0, q_ne_z0,q_np_z0,      Ez0
% initu_func_set = {@initu_func_electrons;@initu_func_ions;0;   0;0;0;    0;0;0};
initu_func_set = {@initu_func_electrons;@initu_func_ions;0;   0;@initq_func_ions_r;0;    0;@initq_func_ions_z;0};

load '../poissonIC350k.mat';
% load '../poissonIC.mat';
UDG_poisson = UDG;      % Load in the poisson as UDG
UDG0 = initu(mesh,initu_func_set,app.arg);      % Change to UDG0 and UH0 for digaso
UDG0(:,[1,2,4,5,7,8],:) = UDG0(:,[1,2,4,5,7,8],:)/ne_star;
UDG0(:,[3,6,9],:) = UDG_poisson;
UH0=inituhat(master,mesh.elcon,UDG0,app.ncu);
[QDG, qq, MiCE] = getq(master, mesh, UDG0, UH0, [], 1);
UDG0(:,app.ncu+1:app.nc,:) = QDG;

% [UDG0,UH0] = getsol(1000);
% UDG0(:,[1,2,4,5,7,8],:) = UDG0(:,[1,2,4,5,7,8],:)/100;
% UH0(1:2,:) = UH0(1:2,:)/100;

% Specific to digaso
app.iterative = 0;
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
app.ncd     = size(mesh.dgnodes,2);
elementtype = mesh.elemtype*ones(mesh.ne,1);
check = 0;

% ------- GMRES Parameters --------- %
% app.restart   = 25;     % 25 is typical 
% app.gmrestol  = 1e-6;   % -6 is a high tolerance
% app.gmresiter = 30;     % 50 is usually the maximum required, for HDG ~10 is typical

app.restart   = 50;     % 25 is typical 
app.gmrestol  = 1e-8;   % -6 is a high tolerance
app.gmresiter = 150;     % 50 is usually the maximum required, for HDG ~10 is typical

% app.restart   = 100;
% app.gmrestol  = 1e-11;
% app.gmresiter = 300;

% ------- Newton Parameters ---------- %
% app.newtoniter = 4;  % def 10 -- if requires >3, then change the timestep
% app.newtontol  = 1e-5; % def 1e-7 -- -5 is fine

app.newtoniter = 6;  % def 10 -- if requires >3, then change the timestep
app.newtontol  = 1e-7; % def 1e-7 -- -5 is fine

% ------- Number of processors ------- %
nproc       = 384;
% nproc       = 4;
app.nfile   = nproc;

% Debug mode
% app.debugmode = 1;

delete *.bin
if nproc>1 % parallel preprocessing
    apppar = digasopre(app,'streamer',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,mesh.bndexpr,[],nproc,0,check,mesh.perm);
    apppar.fileout = 'streamerout';
    return;
else  % serial preprocessing
    appser = digasopre(app,'streamer',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,mesh.bndexpr,[],nproc,0,check,mesh.perm);
    appser.fileout = 'streamerout';
    % Run command:
    % /Users/saustin/Documents/digaso/main/digasoser streamer streamersol
    return;
end
