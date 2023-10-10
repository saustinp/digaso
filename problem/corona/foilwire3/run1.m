setapplicationpath('FM/corona');
 
porder = 2;

DV = 40;
%Va = -80.8;
Vw = DV+Va;
%kappa = 1.355e-3;

scale = 1e3/7*2.59; % geometry scaling factor
Rw=700;
alfa = 160*pi/180;
xc=-209; yc=60;
xw=xc+Rw*cos(alfa);yw=yc+Rw*sin(alfa); rw = 1;

% stabilization parameter
tau = 4;
eps0 = 1;
nu = 1e-3;
t1 = 0;
t2 = 2*pi;
param = {1,nu,eps0,kappa,xw,yw,t1,t2,Wind,tau};

app.uqpk=0;
hybrid = 'hdg';
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fboucorona';
app.fhat = 'fhat';
app.adjoint = 0;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve=1;
app.arg = param;
app.bcm = [8;1;3];
app.bcs = [Vw 0;0 0;Va 0];
app.bcd = [];
app.bcv = []; 

app.denseblock = 0;
app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.fc_q = 1;
app.fc_u = 0;
app.fc_p = 0;

app.np  = 2;
app.nd  = 2;
app.nch = 2;                       % Number of componets of UH
app.nc  = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = 2;

app.appname = 'corona';
app.time = [];
app.dtfc = [];
app.alpha = [];

%mesh = mkmesh_foilwire(porder,400,xw,yw,rw,scale);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% HDG solver
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

etaw = linspace(0,2*pi,50); 
x0 = cos(etaw)*rw+xw;
y0 = sin(etaw)*rw+yw;
nu = 1e-3;
ds = loginc(linspace(0.001,0.12,1000),2);
Ua(:,1,:) = UDG(:,1,:);
Ua(:,2,:) = UDG(:,3,:);
Ua(:,3,:) = UDG(:,5,:);    
[Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh, 0, 1, 0*Ua, Ua, x0, y0, ds, 2);        
figure(2);clf; plot(etaw,Peek1,etaw,R*ones(size(etaw)));

ds = loginc(linspace(0.1,2,500)*7,1);
tetaw = linspace(0,2*pi,400); 
x0 = cos(tetaw)*rw+xw;
y0 = sin(tetaw)*rw+yw;
VDG(:,1,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,5,:)+Wind*mesh.dgnodes(:,4,:);
[x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref); 

return;
