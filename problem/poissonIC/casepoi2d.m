setapplicationpath('FM/axispoi');

% nproc = 1;
fileName = 'poi2d';

porder = 3;
ngrid  = 150;
elemtype = 0;
nodetype = 1;
nstage = 0;
torder = 0;
hybrid = 'hdg';

kappa = 1;  % Diffusivity
l_ref = 1e-4;   % m
mu_ref = 0.05;  % m^2/(V-s)
E_ref = 3e6;    % V/m
e_eps0 = 1.80955e-8;    % Quantity e/epsilon0
phi0 = 18.75e3;     % V
N0 = 5e18; % 1/m^3
z0 = 1e-2; % m
sigma0 = 4e-4;  % m
tau = 1;
app.arg = {kappa, l_ref, mu_ref, E_ref, e_eps0, phi0, N0, z0, sigma0, tau};

app.iterative = 0;
app.getdqdg = 1;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.bcm = [1;2;3;2];
app.bcs = [0;0;0;0];
app.bcd  = [1;1;1;1];  % [nothing,nothing,inlet,outlet,nothing,nothing,airfoil]
app.bcv  = [0;0;0;0];
app.wave = false;
app.tdep = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_u = 1;
app.fc_q = 1;
app.fc_p = 0;
app.nd   = 2;
app.ncu  = 1;               % Number of components of U
app.nch  = app.ncu;                % Number of componets of UH
app.ncq  = app.ncu*app.nd;         % Number of componets of Q
app.nc   = app.ncu+app.ncq;         % Number of componeents of UDG
app.ncp = 0; 
app.time = 0.0;
app.dtfc = 0;
app.alpha = 0;
app.ns   = 1;  
app.adjoint = 0;
app.linear = 1;
app.appname = 'axispoisson';
app.linearSolver = 1;
app.jacobianStep = 0;
app.orderingStep = 0;
app.dt = 0;

% mesh   = mkmesh_square(ngrid,ngrid,porder,0,1,1,elemtype,nodetype);
mesh = mkmesh_rect(100,100,porder,0,[0 125 0 125],0,1);
master = mkmaster(mesh,2*porder);

% preprocessing for c++ code
app.linearproblem = 1;
app.appname='axispoisson';
app.overlappinglevel=1;
app.preconditioner=0;
app.morder = [porder porder];
app.porder = [porder porder];
app.nodetype = nodetype;
app.pgauss = 2*[porder porder];
app.pgaussR = 2*[porder porder];
app.quadtype = [0 0 0];
app.dt = 0;
app.nco = 0;
app.ncd = size(mesh.dgnodes,2);
check = 0;
app.nfile=1;
elementtype = elemtype*ones(mesh.ne,1);
bndexpr=mesh.bndexpr;

[master,mesh,app] = preprocess(master,mesh,app);

UDG0 = initu(mesh,{0;0;0});
UH0 = inituhat(master,mesh.elcon,UDG0,1);

[UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,0*UDG0);
normE = sqrt(UDG(:,2,:).^2 + UDG(:,3,:).^2);

% r = (mesh.dgnodes(:,1,:));
% y = (mesh.dgnodes(:,2,:));
% u = sin(pi*x).*sin(pi*y); % Change this for the axis-symmetry case
% u = exp(-y).*cos(r);
% v = UDG(:,1,:);
% max(abs(u(:)-v(:)))
% return;

% % parallel preprocessing
nproc = 4;
apppar = digasopre(app,'poipar',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,bndexpr,[],nproc,0,check);
apppar.fileout = 'poiparout';

% run this before using plotsol 
system('mpirun -np 4 /Users/saustin/Documents/digaso/main/digasopar poipar poiparout')

[UDGpar,UHpar] = getsolfrombinaryfile(apppar.fileout,apppar.nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
UDGpar = reshape(UDGpar,size(UDG));
normE = sqrt(UDGpar(:,2,:).^2 + UDGpar(:,3,:).^2);

% SERIAL RUN
% nproc = 1;
% app.gmrestol = 1e-12;
% app.newtontol = 1e-12;
% appser = digasopre(app,'poiser',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,bndexpr,[],nproc,0,check);
% appser.fileout = 'poiseroutsol';
% system('/Users/saustin/Documents/digaso/main/digasoser poiser poiserout')
% [UDGser,UHser] = getsolfrombinaryfile(appser.fileout,appser.nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
% UDGser = reshape(UDGser,size(UDG));
% normE = sqrt(UDGser(:,2,:).^2 + UDGser(:,3,:).^2);




% clf;
% figure(); scaplot(mesh,UDGser(:,1,:),[],0,0); axis equal; axis tight;
% figure(); scaplot(mesh,UDG(:,1,:),[],0,0); axis equal; axis tight;

% deformed configuration
% mesh1=mesh;
% mesh1.dgnodes = UDGser(:,1:3,:);
% figure(3); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;



%writeBinaryFile(fileName,mesh,master,app,UDG,UH,[],1);
% app.flag(5) = 2;
% app.factor(5) = 1e-2;
% writeBinaryFile(fileName,mesh,master,app,UDG,UH,[],1);
% 
% 
% [RU,RH,RQ] = hdg_residual(master,app,mesh.dgnodes,mesh.bf,UDG,UH,0*UDG);
% 
% 
% filename = 'solpoi2d';
% fileID = fopen([filename,'.bin'],'r');
% data = fread(fileID,'double');
% fclose(fileID);
% 
% N = master.npv*(mesh.nd+1)*mesh.ne;
% UDG = reshape(data(1:N),[master.npv (mesh.nd+1) mesh.ne]);
% 
% x = (mesh.dgnodes(:,1,:));
% y = (mesh.dgnodes(:,2,:));
% u = sin(pi*x).*sin(pi*y);
% v = UDG(:,1,:);
% max(abs(u(:)-v(:)))
% 
% figure(1); scaplot(mesh,UDG(:,3,:),[],0,1); axis equal; axis tight;
% 
