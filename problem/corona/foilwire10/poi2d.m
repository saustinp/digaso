setapplicationpath('FM/poi');

porder = 2;
elemtype = 0;
nodetype = 1;
hybrid = 'hdg';

kappa = 1;
c = [0,0]; 
tau = 1;

xp=8.5;
the=30*pi/180;
R=5.5;

app1.uqpk=0;
app1.source = 'source';
app1.flux = 'flux';
app1.fbou = 'fbou';
app1.fhat = 'fhat';
app1.localsolve=1;
app1.arg = {kappa,tau};
app1.bcm = [1;1;1];
%app1.bcs = [1;0;0]; % UDGw
app1.bcs = [0;0;1]; % UDGa
app1.bcd = [];
app1.bcv = [];

app1.hybrid = hybrid;
app1.tdep = false;
app1.wave = false;
app1.alag = false;
app1.flg_q = 1;
app1.flg_p = 0;
app1.flg_g = 0;

app1.fc_q = 1;
app1.fc_u = 0;
app1.fc_p = 0;

app1.nd = 2;
app1.nch  = 1;                       % Number of componets of UH
app1.nc   = app1.nch*(app1.nd+1);    % Number of componeents of UDG
app1.ncu = 1;

app1.time = [];
app1.dtfc = [];
app1.alpha = [];

% Rw=700;
% alfa = 160*pi/180;
% xc=-209; yc=60;
%xw=xc+Rw*cos(alfa);yw=yc+Rw*sin(alfa); rw = 1;
%mesh1 = mkmesh1_foilwire(porder,400,xw,yw,rw,scale);
[mesh1,pva,pvw,xw,yw,rw] = mkmesh_foilwire2(porder,3,xp,the,R);
master = mkmaster(mesh1,2*porder);
[master,mesh1] = preprocess(master,mesh1,hybrid);
%mesh1 = periodic(mesh1,{1,'p(:,1)',3,'p(:,1)';4,'p(:,2)',2,'p(:,2)'});

UDGp = initu(mesh1,{0;0;0});
UHp=inituhat(master,mesh1.elcon,UDGp,1);

% HDG solver
[UDGp,UHp] = hdg_solve(master,mesh1,app1,UDGp,UHp,0*UDGp);

% mesh11 = mkmesh1_wireplate(porder+1,0);
% master1 = mkmaster(mesh11,2*(porder+1));
% [master1,mesh11] = preprocess(master1,mesh11,hybrid);
% UDGpstar = postprocessnd(master,mesh1,master1,mesh11,UDGp);

figure(1); clf; scaplot(mesh1,UDGp(:,1,:),[],1,0); 
axis equal; axis tight; colormap jet;

% figure(2); clf; scaplot(mesh1,sqrt(UDGp(:,2,:).^2+UDGp(:,3,:).^2),[],1,0); 
% axis equal; axis tight; colormap jet;

% x = mesh1.p(:,1);
% y = mesh1.p(:,2);
% [Vx, Vy, Psix, Psiy, Psi, Psiax, Psiay, Psia, Vax, Vay] = getvelocity(x, y, xw, yw, rw, scale);
% t = mesh1.t'; t(4,:) = 1; p = mesh1.p';
% figure(1); clf; pdeplot(p,[],t,'XYData',Psi(:),'Contour','on','Levels',100);  axis tight; axis on; colormap jet;

return;

porder=2;
mesh1 = mkmesh1_foilwire(porder,400,xw,yw,rw,scale);

Vw = 10;
Va = -10;
[Ex, Ey, U] = getfield(mesh1.dgnodes(:,1,:), mesh1.dgnodes(:,2,:), Va, Vw, xw, yw, rw, [], [], [], scale);
UDG = UDGa*Va + UDGw*Vw;
figure(1); clf; scaplot(mesh1,U,[-10 10],1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDG(:,1,:),[-10 10],1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh1,Ex,[-10 10]/100,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDG(:,2,:),[-10 10]/100,1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh1,Ey,[-10 10]/40,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDG(:,3,:),[-10 10]/40,1,0); axis equal; axis tight; colormap jet;

[Vx, Vy, Psix, Psiy, Psi, Psiax, Psiay, Psia, Vax, Vay] = getvelocity(mesh1.dgnodes(:,1,:), mesh1.dgnodes(:,2,:), xw, yw, rw, scale);
x=mesh1.dgnodes(:,1,:);
[~,imin] = min(x(:));
Vx = Vx/Vx(imin);
mesh1.dgnodes(:,3,:) = Vx;
mesh1.dgnodes(:,4,:) = Vy;

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
VDG(:,1,:) = UDG(:,2,:)+Wind*mesh1.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,3,:)+Wind*mesh1.dgnodes(:,4,:);
[x,y,s,Ex,Ey,in,out] = fieldlines(mesh1, VDG, x0, y0, ds, nref); 
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
VDG(:,1,:) = UDG(:,2,:)+Wind*mesh1.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,3,:)+Wind*mesh1.dgnodes(:,4,:);
[x,y,s,Ex,Ey,in,out] = fieldlines(mesh1, VDG, x0, y0, ds, nref); 


save poi.mat UDGpoi1 UDGpoi2;

for i = 1:3
    figure(i); clf; scaplot(mesh1,UDGpoi1(:,i,:),[],1,0); 
    axis equal; axis tight; colormap jet;
end

[Ex, Ey, U] = getfield(mesh1.dgnodes(:,1,:), mesh1.dgnodes(:,2,:), -10, 10, xw, yw, rw, [], [], [], scale);
figure(1); clf; scaplot(mesh1,U,[-10 10],1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDGp(:,1,:),[-10 10],1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh1,Ex,[-10 10]/400,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDGp(:,2,:),[-10 10]/400,1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh1,Ey,[-10 10]/40,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDGp(:,3,:),[-10 10]/40,1,0); axis equal; axis tight; colormap jet;

[Ex1, Ey1, U1] = getfield(mesh1.dgnodes(:,1,:), mesh1.dgnodes(:,2,:), 0, 1, xw, yw, rw, [], [], [], scale);
[Ex2, Ey2, U2] = getfield(mesh1.dgnodes(:,1,:), mesh1.dgnodes(:,2,:), 1, 0, xw, yw, rw, [], [], [], scale);

figure(1); clf; scaplot(mesh1,U1,[0 1],1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDGpoi1(:,1,:),[0 1],1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh1,Ex1,[-10 10]/400,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDGpoi1(:,2,:),[-10 10]/400,1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh1,Ey1,[-10 10]/40,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDGpoi1(:,3,:),[-10 10]/40,1,0); axis equal; axis tight; colormap jet;

figure(1); clf; scaplot(mesh1,U2,[0 1],1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDGpoi2(:,1,:),[0 1],1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh1,Ex2,[-10 10]/4000,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDGpoi2(:,2,:),[-10 10]/4000,1,0); axis equal; axis tight; colormap jet;
figure(1); clf; scaplot(mesh1,Ey2,[-10 10]/4000,1,0); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh1,UDGpoi2(:,3,:),[-10 10]/4000,1,0); axis equal; axis tight; colormap jet;


nn = [2 4 8 16 32]+1;
for porder = 1:4
    for ii=1:length(nn)
        ngrid = nn(ii);
        poi2d;
        e = calerror(UDGp,mesh1,master,@exactsol1);          
        erru(ii,porder) = e(1); errq(ii,porder) = sqrt(e(2)^2+e(3)^2);
        errs(ii,porder) = calerror(UDGpstar,mesh11,master1,@exactsol1);
    end
end

cu=log(erru(1:end-1,:)./erru(2:end,:))/log(2);
cs=log(errs(1:end-1,:)./errs(2:end,:))/log(2);
cq=log(errq(1:end-1,:)./errq(2:end,:))/log(2);
a=[erru [0 0 0 0; cu]; errq [0 0 0 0; cq]; errs [0 0 0 0; cs]];
a=a(:,[1 5 2 6 3 7 4 8]);
a=[[nn-1 nn-1 nn-1]' a];


figure(1); clf; scaplot3(mesh1,UDGp(:,1,:),[-0.2 1.2],2); 
%scaplot(mesh1,UDGp(:,1,:),[0 1],1,1,1)
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

figure(1); clf; scaplot3(mesh11,UDGpstar(:,1,:),[-0.2 1.2],2); 
%scaplot(mesh1,UDGp(:,1,:),[0 1],1,1,1)
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

[RU,RH,RQ] = hdg_residual(master,app,mesh1,UDGp,UHp,0*UDGp);

syms x y 
u = sin(0.5*pi*x).*sin(0.5*pi*y);
ux = diff(u,'x');
uxx = diff(ux,'x');
uy = diff(u,'y');
uyy = diff(uy,'y');
f = -simplify(uxx+uyy);

% app1.ubou = 'ldgubou';
% UHp0 = ldg_uhat(mesh1,master,app,UDGp(:,1,:));

% app1.source = 'source';
% app1.flux = 'flux';
% app1.fbou = 'ldgfbou';
% app1.fhat = 'ldgfhat';
%[RU0,UDGp0,UHp0] = ldg_residual(master,app,mesh1,UDGp(:,1,:),0*UDGp(:,1,:));

% app1.source = 'source';
% app1.flux = 'flux';
% app1.fbou = 'fbou';
% app1.fhat = 'fhat';
% [RU,RH,RQ] = hdg_residual(master,app,mesh1,UDGp,UHp,0*UDGp);

% 
app1.source = 'source';
app1.flux = 'flux';
app1.ubou = 'ldgubou';
app1.fbou = 'ldgfbou';
app1.fhat = 'ldgfhat';
tol = 1e-8;
dt = 1e-1*ones(1000,1);
[UDGpA,UHpA,normR] = ardm(master,app,mesh1,0*UDGp(:,1,:),dt,tol);

figure(1); scaplot(mesh1,UDGpA(:,1,:),[],0,1); axis equal; axis tight; colormap jet;

% [un,normR] = restartedardm(master,app,mesh1,UDGp,UHp,0*SH,dt,tol);


% % HDG postprocessing 
% mesh11 = mkmesh1_square(ngrid,ngrid,porder+1,0,1,1,elemtype,nodetype);
% master1 = mkmaster(mesh11,2*(porder+1));
% [master1,mesh11] = preprocess(master1,mesh11,hybrid);
% UDGpstar = postprocessnd(master,mesh1,master1,mesh11,UDGp);
% 
% figure(2); scaplot(mesh11,UDGpstar(:,1,:),[],0,1); axis equal; axis tight;
% %figure(2); scaplot(mesh1,u(:,1,:),[],0,1); axis equal; axis tight;
% 
% x = (mesh11.dgnodes(:,1,:));
% y = (mesh11.dgnodes(:,2,:));
% u = sin(pi*x).*sin(pi*y);
% v = UDGpstar(:,1,:);
% max(abs(u(:)-v(:)))
% 
