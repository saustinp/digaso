setapplicationpath('FM/poi');

porder = 1;
ngrid  = 3;
elemtype = 0;
nodetype = 1;
hybrid = 'hdg';

kappa = 1;
c = [0,0]; 
tau = 1;

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.localsolve=1;
app.arg = {kappa,tau};
app.bcm = [5];
app.bcs = [0];
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

mesh = mkmesh_Lshape(porder,ngrid,elemtype,nodetype);
master = mkmaster(mesh,2*(porder+1));
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

% HDG solver
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,0*UDG);

mesh1 = mkmesh_Lshape(porder+1,ngrid,elemtype,nodetype);
master1 = mkmaster(mesh1,2*(porder+2));
[master1,mesh1] = preprocess(master1,mesh1,hybrid);
UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);

%[UDG,UH]=hdg_solve(master,mesh,app,UDG,UH,SH)
% figure(1); clf; scaplot(mesh,UDG(:,1,:),[],0,1); axis equal; axis tight; colormap jet;
% 
% x = (mesh.dgnodes(:,1,:));
% y = (mesh.dgnodes(:,2,:));
% r = sqrt(x.^2 + y.^2);
% t = atan(y./x);
% i = (x==0 & y==0);
% t(i) = 0;
% i = (x<0);
% t(i) = t(i)+pi;
% i = (x==0 & y<0);
% t(i) = 3*pi/2;
% i = (x>0 & y<0);
% t(i) = t(i)+2*pi;
% u = sin((2*t)/3).*((r).^(2/3));

% v = UDG(:,1,:);
% max(abs(u(:)-v(:)))
return;


nn = [2 4 8 16 32]+1;
for porder = 1:3
    for ii=1:length(nn)
        ngrid = nn(ii);
        poi2d;
        e = calerror(UDG,mesh,master,@exactsol1);          
        erru(ii,porder) = e(1); errq(ii,porder) = sqrt(e(2)^2+e(3)^2);
        errs(ii,porder) = calerror(UDGstar,mesh1,master1,@exactsol1);
    end
end

cu=log(erru(1:end-1,:)./erru(2:end,:))/log(2);
cs=log(errs(1:end-1,:)./errs(2:end,:))/log(2);
cq=log(errq(1:end-1,:)./errq(2:end,:))/log(2);
a=[erru [0 0 0 0; cu]; errq [0 0 0 0; cq]; errs [0 0 0 0; cs]];
a=a(:,[1 5 2 6 3 7]);
a=[[nn-1 nn-1 nn-1]' a];


figure(1); clf; %scaplot3(mesh,UDG(:,1,:),[-0.2 1.2],2); 
scaplot(mesh,UDG(:,1,:),[-0.2 1.4],1,1)
axis equal; axis tight; colormap jet;
axis normal; 
set(gca,'FontSize',16);
set(gca,'xtick',[-1:0.2:1]);
set(gca,'ytick',[-1:0.2:1]);
%set(gca,'ztick',[-0.2:0.2:142]);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
box on;


figure(1); clf; %scaplot3(mesh,UDG(:,1,:),[-0.2 1.2],2); 
scaplot(mesh1,UDGstar(:,1,:),[-0.2 1.4],1,1)
axis equal; axis tight; colormap jet;
axis normal; 
set(gca,'FontSize',16);
set(gca,'xtick',[-1:0.2:1]);
set(gca,'ytick',[-1:0.2:1]);
%set(gca,'ztick',[-0.2:0.2:142]);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
box on;

figure(1); clf; scaplot3(mesh1,UDGstar(:,1,:),[-0.2 1.2],2); 
%scaplot(mesh,UDG(:,1,:),[0 1],1,1,1)
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

[RU,RH,RQ] = hdg_residual(master,app,mesh,UDG,UH,0*UDG);

syms x y 
u = sin(0.5*pi*x).*sin(0.5*pi*y);
ux = diff(u,'x');
uxx = diff(ux,'x');
uy = diff(u,'y');
uyy = diff(uy,'y');
f = -simplify(uxx+uyy);

% app.ubou = 'ldgubou';
% UH0 = ldg_uhat(mesh,master,app,UDG(:,1,:));

% app.source = 'source';
% app.flux = 'flux';
% app.fbou = 'ldgfbou';
% app.fhat = 'ldgfhat';
%[RU0,UDG0,UH0] = ldg_residual(master,app,mesh,UDG(:,1,:),0*UDG(:,1,:));

% app.source = 'source';
% app.flux = 'flux';
% app.fbou = 'fbou';
% app.fhat = 'fhat';
% [RU,RH,RQ] = hdg_residual(master,app,mesh,UDG,UH,0*UDG);

% 
app.source = 'source';
app.flux = 'flux';
app.ubou = 'ldgubou';
app.fbou = 'ldgfbou';
app.fhat = 'ldgfhat';
tol = 1e-8;
dt = 1e-1*ones(1000,1);
[UDGA,UHA,normR] = ardm(master,app,mesh,0*UDG(:,1,:),dt,tol);

figure(1); scaplot(mesh,UDGA(:,1,:),[],0,1); axis equal; axis tight; colormap jet;

% [un,normR] = restartedardm(master,app,mesh,UDG,UH,0*SH,dt,tol);


% % HDG postprocessing 
% mesh1 = mkmesh_square(ngrid,ngrid,porder+1,0,1,1,elemtype,nodetype);
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);
% 
% figure(2); scaplot(mesh1,UDGstar(:,1,:),[],0,1); axis equal; axis tight;
% %figure(2); scaplot(mesh,u(:,1,:),[],0,1); axis equal; axis tight;
% 
% x = (mesh1.dgnodes(:,1,:));
% y = (mesh1.dgnodes(:,2,:));
% u = sin(pi*x).*sin(pi*y);
% v = UDGstar(:,1,:);
% max(abs(u(:)-v(:)))
% 
