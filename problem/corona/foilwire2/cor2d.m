setapplicationpath('FM/corona');

% See the paper Neimarlija2009.pdf and Medlin1998b.pdf
 
porder = 2;


DV = 40;
Va = -68;
Vw = DV+Va;

scale = 1e3/7*2.59; % geometry scaling factor
Rw=900;
alfa = 5*pi/180;
xc=-209; yc=60;
xw=xc+Rw*cos(alfa);yw=yc+Rw*sin(alfa); rw = 1;

% stabilization parameter
tau = 4;
eps0 = 1;
kappa = 1.25e-2/9;
nu = 1e-3;
% t1 = 3.4;
% t2 = 5.4;
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

% [Vx, Vy, Psix, Psiy, Psi, Psiax, Psiay, Psia, Vax, Vay] = getvelocity(mesh.dgnodes(:,1,:), mesh.dgnodes(:,2,:), xw, yw, rw, scale);
% x=mesh.dgnodes(:,1,:);
% [~,imin] = min(x(:));
% Vx = Vx/Vx(imin);
% mesh.dgnodes(:,3,:) = Vx;
% mesh.dgnodes(:,4,:) = Vy;

% UDG = initu(mesh,{1;0;0;0;0;0});
% Ua = Vw*UDGw + Va*UDGa;
% UDG(:,1,:) = Ua(:,1,:);
% UDG(:,3,:) = Ua(:,2,:);
% UDG(:,5,:) = Ua(:,3,:);
% UH = inituhat(master,mesh.elcon,UDG,app.ncu);

% HDG solver
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

for i = 1:length(jj)
    figure(i); clf; scaplot(mesh,UDG(:,jj(i),:),clm{jj(i)},2); 
    axis equal; axis([-800 1600 -800 1200]); axis off; colormap jet; colorbar('FontSize',15);    
    %set(gca, 'LooseInset', get(gca, 'TightInset'));
%     fn = ['numsol0' num2str(i)];
%     print('-dpng',fn);
end

etaw = linspace(0,2*pi,50); 
x0 = cos(etaw)*rw+xw;
y0 = sin(etaw)*rw+yw;
nu = 1e-3;
ds = loginc(linspace(0.002,0.1,500),2);
Ua(:,1,:) = UDG(:,1,:);
Ua(:,2,:) = UDG(:,3,:);
Ua(:,3,:) = UDG(:,5,:);    
[Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh, 0, 1, 0*Ua, Ua, x0, y0, ds, 2);        
figure(2);clf; plot(etaw,Peek1,etaw,R*ones(size(etaw)));

ds = loginc(linspace(0.1,1,500)*4,1);
tetaw = linspace(0,2*pi,200); 
x0 = cos(tetaw)*rw+xw;
y0 = sin(tetaw)*rw+yw;
VDG(:,1,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,5,:)+Wind*mesh.dgnodes(:,4,:);
[x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref); 

return;



% for i = 1:6
%     figure(i); clf; scaplot(mesh,UDG(:,i,:),[],2); 
%     axis equal; axis tight; axis off; colormap jet;
%     colorbar('FontSize',15);
%     hold on; 
%     plot(xw,yw,'o');
%     %set(gca, 'LooseInset', get(gca, 'TightInset'));
%     fn = ['numsol' num2str(i)];
%     print('-dpng',fn);
% end

nu=nu/100;
param = {1,nu,eps0,kappa,xw,yw,t1,t2,tau};
app.arg = param;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

clm{1} = [];
clm{2} = [];
clm{3} = [-1 1]/5;
clm{5} = [-1 1]/5;
jj = [1 2 3 5];
for i = 1:length(jj)
    figure(i); clf; scaplot(mesh,UDG(:,jj(i),:),clm{jj(i)},2); 
    axis equal; axis([-800 1600 -800 1200]); axis off; colormap jet; colorbar('FontSize',15);    
    %set(gca, 'LooseInset', get(gca, 'TightInset'));
%     fn = ['numsol0' num2str(i)];
%     print('-dpng',fn);
end
figure(5); clf; scaplot(mesh,mesh.dgnodes(:,3,:),[-1 1]/5,2); 
axis equal; axis([-800 1000 -400 800]); axis off; colormap jet; colorbar('FontSize',15);    
figure(6); clf; scaplot(mesh,mesh.dgnodes(:,4,:),[-1 1]/5,2); 
axis equal; axis([-800 1000 -400 800]); axis off; colormap jet; colorbar('FontSize',15);    
figure(7); clf; scaplot(mesh,UDG(:,3,:)+mesh.dgnodes(:,3,:),[-1 1]/5,2); 
axis equal; axis([-800 1000 -400 800]); axis off; colormap jet; colorbar('FontSize',15);    
figure(8); clf; scaplot(mesh,UDG(:,5,:)+mesh.dgnodes(:,4,:),[-1 1]/5,2); 
axis equal; axis([-800 1000 -400 800]); axis off; colormap jet; colorbar('FontSize',15);    


ds = linspace(0.1,1,1000)*2;
VDG(:,1,:) = UDG(:,3,:)+mesh.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,5,:)+mesh.dgnodes(:,4,:);
tetaw = linspace(t1,t2,50); 
x1 = cos(tetaw)*rw+xw;
y1 = sin(tetaw)*rw+yw;
[x,y,s,Ex,Ey] = fieldlines(mesh, VDG, x1, y1, ds, 1);

param = {1,1e-1/100,eps0,kappa,xw,yw,t1,t2,tau};
app.arg = param;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

Ua(:,1,:) = UDG(:,1,:);
Ua(:,2,:) = UDG(:,3,:);
Ua(:,3,:) = UDG(:,5,:);
%ds = linspace(0.01,0.1,5000)/50;
ds = logdec(linspace(1e-3,0.1,200),2);
[Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh, 0, 1, 0*Ua, Ua, x0, y0, ds, 2);

[~,ind] = max(Peek1);
figure(1);clf;plot(s1(:,ind),sqrt(Ex1(:,ind).^2+Ey1(:,ind).^2));

[~,ind] = sort(abs(Peek1-R));
figure(4);clf;
plot(etaw,Peek1,etaw,R*ones(size(etaw)));
% hold on;
% plot([etaw(ind(1)) etaw(ind(1))], [min(Peek) R], '-k');
% plot([etaw(ind(2)) etaw(ind(2))], [min(Peek) R], '-k');
xlabel('\theta','FontSize',18);
ylabel('Peek Integral','FontSize',18);
title('\Delta V = 38, V_a = -14, V_w = 24','FontSize',18);
set(gca,'FontSize',16);
axis tight;

ds = linspace(0.01,0.1,2000)/20;
[Peek,x,y,s,Ex,Ey] = PeekIntegral(mesh, Vw, Va, UDGw, UDGa, x0, y0, ds, 1);
save boundata2.mat etaw Peek xw yw R

[~,ind] = sort(abs(Peek-R));
figure(4);clf;
plot(etaw,Peek,etaw,R*ones(size(etaw)));
hold on;
plot([etaw(ind(1)) etaw(ind(1))], [min(Peek) R], '-k');
plot([etaw(ind(2)) etaw(ind(2))], [min(Peek) R], '-k');
xlabel('\theta','FontSize',18);
ylabel('Peek Integral','FontSize',18);
title('\Delta V = 34, V_a = -9, V_w = 25','FontSize',18);
set(gca,'FontSize',16);
axis tight;


t1 = 3.2;
t2 = 4.2;
tm = 0.5*(t1+t2);
td = t2-t1;
tc = (log(1e-3))/(abs(t1-tm).^3);
tt = linspace(0,2*pi,1000);
ft = exp(tc*abs(tt-tm).^3);
ft(tt>t2)=0;
ft(tt<t1)=0;
figure(1);clf;plot(tt,ft);









