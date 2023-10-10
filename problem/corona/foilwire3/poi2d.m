setapplicationpath('FM/poi');

porder = 2;
elemtype = 0;
nodetype = 1;
hybrid = 'hdg';

kappa = 1;
c = [0,0]; 
tau = 1;

app.uqpk=0;
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.localsolve=1;
app.arg = {kappa,tau};
app.bcm = [1;1;1];
%app.bcs = [1;0;0]; % UDGw
app.bcs = [0;0;1]; % UDGa
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
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = 1;

app.time = [];
app.dtfc = [];
app.alpha = [];

scale = 1e3/7*2.59; % geometry scaling factor
Rw=700;
alfa = 160*pi/180;
xc=-209; yc=60;
xw=xc+Rw*cos(alfa);yw=yc+Rw*sin(alfa); rw = 1;
mesh = mkmesh_foilwire(porder,400,xw,yw,rw,scale);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
%mesh = periodic(mesh,{1,'p(:,1)',3,'p(:,1)';4,'p(:,2)',2,'p(:,2)'});

UDGp = initu(mesh,{0;0;0});
UHp=inituhat(master,mesh.elcon,UDGp,1);

% HDG solver
[UDGp,UHp] = hdg_solve(master,mesh,app,UDGp,UHp,0*UDGp);

% mesh1 = mkmesh_wireplate(porder+1,0);
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% UDGpstar = postprocessnd(master,mesh,master1,mesh1,UDGp);

figure(1); clf; scaplot(mesh,UDGp(:,1,:),[],1,0); 
axis equal; axis tight; colormap jet;

% figure(2); clf; scaplot(mesh,sqrt(UDGp(:,2,:).^2+UDGp(:,3,:).^2),[],1,0); 
% axis equal; axis tight; colormap jet;

% x = mesh.p(:,1);
% y = mesh.p(:,2);
% [Vx, Vy, Psix, Psiy, Psi, Psiax, Psiay, Psia, Vax, Vay] = getvelocity(x, y, xw, yw, rw, scale);
% t = mesh.t'; t(4,:) = 1; p = mesh.p';
% figure(1); clf; pdeplot(p,[],t,'XYData',Psi(:),'Contour','on','Levels',100);  axis tight; axis on; colormap jet;

return;

porder=2;
mesh = mkmesh_foilwire(porder,400,xw,yw,rw,scale);

Vw = 10;
Va = -10;
[Ex, Ey, U] = getfield(mesh.dgnodes(:,1,:), mesh.dgnodes(:,2,:), Va, Vw, xw, yw, rw, [], [], [], scale);
UDG = UDGa*Va + UDGw*Vw;
figure(1); clf; scaplot(mesh,U,[-10 10],1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDG(:,1,:),[-10 10],1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh,Ex,[-10 10]/100,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDG(:,2,:),[-10 10]/100,1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh,Ey,[-10 10]/40,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDG(:,3,:),[-10 10]/40,1,0); axis equal; axis tight; colormap jet;

[Vx, Vy, Psix, Psiy, Psi, Psiax, Psiay, Psia, Vax, Vay] = getvelocity(mesh.dgnodes(:,1,:), mesh.dgnodes(:,2,:), xw, yw, rw, scale);
x=mesh.dgnodes(:,1,:);
[~,imin] = min(x(:));
Vx = Vx/Vx(imin);
mesh.dgnodes(:,3,:) = Vx;
mesh.dgnodes(:,4,:) = Vy;

nref=2;
DV=35;
Va=-20;
Vw=Va+DV;
Wind = 0.0;
etaw = linspace(0,2*pi,160);
x0 = cos(etaw)*rw+xw;
y0 = sin(etaw)*rw+yw;
ds = linspace(1e-1,1,1000)*2;
UDG = UDGa*Va + UDGw*Vw;
VDG(:,1,:) = UDG(:,2,:)+Wind*mesh.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,4,:);
[x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref); 
[x,y,s,Ex,Ey] = streamlines(Va, DV, Wind, xw, yw, rw, [], [], [], x0, y0, ds, scale);

figure(1); axis([-500 700 -300 800]);
figure(2); axis([-500 700 -200 800]);


DV=35;
Va=-20;
Vw=Va+DV;
Wind = 0.05;
etaw = linspace(3.4,5.4,60);
x0 = cos(etaw)*rw+xw;
y0 = sin(etaw)*rw+yw;
ds = linspace(1e-1,1,1000)*2;
UDG = UDGa*Va + UDGw*Vw;
VDG(:,1,:) = UDG(:,2,:)+Wind*mesh.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,4,:);
[x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref); 


save poi.mat UDGpoi1 UDGpoi2;

for i = 1:3
    figure(i); clf; scaplot(mesh,UDGpoi1(:,i,:),[],1,0); 
    axis equal; axis tight; colormap jet;
end

[Ex, Ey, U] = getfield(mesh.dgnodes(:,1,:), mesh.dgnodes(:,2,:), -10, 10, xw, yw, rw, [], [], [], scale);
figure(1); clf; scaplot(mesh,U,[-10 10],1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDGp(:,1,:),[-10 10],1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh,Ex,[-10 10]/400,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDGp(:,2,:),[-10 10]/400,1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh,Ey,[-10 10]/40,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDGp(:,3,:),[-10 10]/40,1,0); axis equal; axis tight; colormap jet;

[Ex1, Ey1, U1] = getfield(mesh.dgnodes(:,1,:), mesh.dgnodes(:,2,:), 0, 1, xw, yw, rw, [], [], [], scale);
[Ex2, Ey2, U2] = getfield(mesh.dgnodes(:,1,:), mesh.dgnodes(:,2,:), 1, 0, xw, yw, rw, [], [], [], scale);

figure(1); clf; scaplot(mesh,U1,[0 1],1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDGpoi1(:,1,:),[0 1],1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh,Ex1,[-10 10]/400,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDGpoi1(:,2,:),[-10 10]/400,1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh,Ey1,[-10 10]/40,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDGpoi1(:,3,:),[-10 10]/40,1,0); axis equal; axis tight; colormap jet;

figure(1); clf; scaplot(mesh,U2,[0 1],1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDGpoi2(:,1,:),[0 1],1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh,Ex2,[-10 10]/4000,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDGpoi2(:,2,:),[-10 10]/4000,1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh,Ey2,[-10 10]/4000,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,UDGpoi2(:,3,:),[-10 10]/4000,1,0); axis equal; axis tight; colormap jet;


nn = [2 4 8 16 32]+1;
for porder = 1:4
    for ii=1:length(nn)
        ngrid = nn(ii);
        poi2d;
        e = calerror(UDGp,mesh,master,@exactsol1);          
        erru(ii,porder) = e(1); errq(ii,porder) = sqrt(e(2)^2+e(3)^2);
        errs(ii,porder) = calerror(UDGpstar,mesh1,master1,@exactsol1);
    end
end

cu=log(erru(1:end-1,:)./erru(2:end,:))/log(2);
cs=log(errs(1:end-1,:)./errs(2:end,:))/log(2);
cq=log(errq(1:end-1,:)./errq(2:end,:))/log(2);
a=[erru [0 0 0 0; cu]; errq [0 0 0 0; cq]; errs [0 0 0 0; cs]];
a=a(:,[1 5 2 6 3 7 4 8]);
a=[[nn-1 nn-1 nn-1]' a];


figure(1); clf; scaplot3(mesh,UDGp(:,1,:),[-0.2 1.2],2); 
%scaplot(mesh,UDGp(:,1,:),[0 1],1,1,1)
axis equal; axis tight; colormap jet;
axis([0 1 0 1 -0.2 1.2]);
axis normal; colorbar off;
set(gca,'FontSize',16);
set(gca,'xtick',[0:0.2:1]);
set(gca,'ytick',[0:0.2:1]);
set(gca,'ztick',[-0.2:0.2:1.2]);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
box on;

figure(1); clf; scaplot3(mesh1,UDGpstar(:,1,:),[-0.2 1.2],2); 
%scaplot(mesh,UDGp(:,1,:),[0 1],1,1,1)
axis equal; axis tight; colormap jet;
axis([0 1 0 1 -0.2 1.2]);
axis normal; colorbar off;
set(gca,'FontSize',16);
set(gca,'xtick',[0:0.2:1]);
set(gca,'ytick',[0:0.2:1]);
set(gca,'ztick',[-0.2:0.2:1.2]);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
box on;

[RU,RH,RQ] = hdg_residual(master,app,mesh,UDGp,UHp,0*UDGp);

syms x y 
u = sin(0.5*pi*x).*sin(0.5*pi*y);
ux = diff(u,'x');
uxx = diff(ux,'x');
uy = diff(u,'y');
uyy = diff(uy,'y');
f = -simplify(uxx+uyy);

% app.ubou = 'ldgubou';
% UHp0 = ldg_uhat(mesh,master,app,UDGp(:,1,:));

% app.source = 'source';
% app.flux = 'flux';
% app.fbou = 'ldgfbou';
% app.fhat = 'ldgfhat';
%[RU0,UDGp0,UHp0] = ldg_residual(master,app,mesh,UDGp(:,1,:),0*UDGp(:,1,:));

% app.source = 'source';
% app.flux = 'flux';
% app.fbou = 'fbou';
% app.fhat = 'fhat';
% [RU,RH,RQ] = hdg_residual(master,app,mesh,UDGp,UHp,0*UDGp);

% 
app.source = 'source';
app.flux = 'flux';
app.ubou = 'ldgubou';
app.fbou = 'ldgfbou';
app.fhat = 'ldgfhat';
tol = 1e-8;
dt = 1e-1*ones(1000,1);
[UDGpA,UHpA,normR] = ardm(master,app,mesh,0*UDGp(:,1,:),dt,tol);

figure(1); scaplot(mesh,UDGpA(:,1,:),[],0,1); axis equal; axis tight; colormap jet;

% [un,normR] = restartedardm(master,app,mesh,UDGp,UHp,0*SH,dt,tol);


% % HDG postprocessing 
% mesh1 = mkmesh_square(ngrid,ngrid,porder+1,0,1,1,elemtype,nodetype);
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% UDGpstar = postprocessnd(master,mesh,master1,mesh1,UDGp);
% 
% figure(2); scaplot(mesh1,UDGpstar(:,1,:),[],0,1); axis equal; axis tight;
% %figure(2); scaplot(mesh,u(:,1,:),[],0,1); axis equal; axis tight;
% 
% x = (mesh1.dgnodes(:,1,:));
% y = (mesh1.dgnodes(:,2,:));
% u = sin(pi*x).*sin(pi*y);
% v = UDGpstar(:,1,:);
% max(abs(u(:)-v(:)))
% 
