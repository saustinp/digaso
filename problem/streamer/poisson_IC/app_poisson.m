porder = 2;
hybrid = 'hdg';
kappa = 1;
tau = 1;
clear app;
app.axisymmetry = 1;

app.arg = phys_param();     % Physics param loaded in a separate script
app.arg{end+1} = tau;

app.source = 'source_poisson';
app.flux = 'flux_poisson';
app.fbou = 'fbou_poisson';
app.fhat = 'fhat_poisson';
app.localsolve=1;
app.bcm = [1;3;1;3];

l_ref = app.arg{2};
E_ref = app.arg{4};
phi0 = app.arg{6};
phi0_tilde = phi0/(E_ref*l_ref);

app.bcs = [0;0;phi0_tilde;0];
% app.bcs = [0;0;0;0];
app.bcd = [];
app.bcv = [];

app.hybrid = hybrid;
app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.fc_q = 1;
app.fc_u = 0;
app.fc_p = 0;

app.nd = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of components of UDG: for each equation, one solution variable and ND number of gradient variables
app.ncu = 1;

app.time = [];
app.dtfc = [];
app.alpha = [];

% Initializing data structures
% mesh = mkmesh_rect(41,81,porder,0,[0 125 0 125],0,1);
mesh = mkmesh_streamer_gmsh(porder, "streamer_16k_fixed.msh");
% mesh = mkmesh_streamer_gmsh(porder, "streamer_57k.msh");
% mesh = mkmesh_streamer_gmsh(porder, "streamer_89k.msh");
% mesh = mkmesh_streamer_gmsh(porder, "streamer_114k.msh");

master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,0*UDG);

save '../poissonIC.mat' UDG

normE = sqrt(UDG(:,2,:).^2 + UDG(:,3,:).^2);
figure(); scaplot(mesh,normE,[],0,0); axis equal; axis tight; colormap jet;
