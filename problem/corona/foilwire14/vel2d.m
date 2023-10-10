setapplicationpath('FM/poi');

porder = 2;
elemtype = 0;
nodetype = 1;
hybrid = 'hdg';

kappa = 1;
c = [0,0]; 
tau = 1;

% xp=8.5;
% the=30*pi/180;
% R=5.5;

app1.uqpk=0;
app1.source = 'source';
app1.flux = 'flux';
app1.fbou = 'fbou';
app1.fhat = 'fhat';
app1.localsolve=1;
app1.arg = {kappa,tau};
app1.bcm = [3;4;3];
app1.bcs = [0;0;0];
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

[mesh1,pva,pvw,xw,yw,rw,pvc] = mkmesh_foilwire2(porder,3,xp,the,R);
master = mkmaster(mesh1,2*porder);
[master,mesh1] = preprocess(master,mesh1,hybrid);

UDGp = initu(mesh1,{0;0;0});
UHp=inituhat(master,mesh1.elcon,UDGp,1);

% HDG solver
[UDGp,UHp] = hdg_solve(master,mesh1,app1,UDGp,UHp,0*UDGp);

figure(1); clf; scaplot(mesh1,UDGp(:,1,:),[],1,0); 
axis equal; axis tight; colormap jet;
hold on; 
plot(pva(:,1),pva(:,2),'-k','LineWidth',1.5);
plot(pvw(:,1),pvw(:,2),'-k','LineWidth',1.5);

figure(2); clf; scaplot(mesh1,-UDGp(:,2,:),[0 1.5],1);
axis equal; axis tight; colormap jet;
hold on; 
plot(pva(:,1),pva(:,2),'-k','LineWidth',1.5);
plot(pvw(:,1),pvw(:,2),'-k','LineWidth',1.5);

figure(3); clf; scaplot(mesh1,-UDGp(:,3,:),[-.25 .25],1);
axis equal; axis tight; colormap jet;
hold on; 
plot(pva(:,1),pva(:,2),'-ko','LineWidth',1.5);
plot(pvw(:,1),pvw(:,2),'-ko','LineWidth',1.5);

mesh1.dgnodes(:,3,:) = -UDGp(:,2,:);
mesh1.dgnodes(:,4,:) = -UDGp(:,3,:);

return;

[Psi] = fieldatx(mesh1,UDGp(:,1,:),mesh1.p,2);
t = mesh1.t'; t(4,:) = 1; p = mesh1.p'; 
figure(4); clf; pdeplot(p,[],t,'XYData',Psi(:),'Contour','on','Levels',100);  
axis tight; axis on; colormap jet;
hold on; plot(pvw(:,1),pvw(:,2),'-ko','LineWidth',1.5);

U = fieldatx(mesh1,mesh1.dgnodes(:,3,:),mesh1.p,2);
V = fieldatx(mesh1,mesh1.dgnodes(:,4,:),mesh1.p,2);
t = mesh1.t'; t(4,:) = 1; p = mesh1.p'; 
figure(5); clf; hold on;
pdeplot(p,[],t,'FlowData',[U(:) V(:)],'FlowStyle','off','mesh','on');  
hold on;
pdeplot(p,[],t,'FlowData',[U(:) V(:)],'FlowStyle','arrow','mesh','on');  
plot(pva(:,1),pva(:,2),'-k','LineWidth',1.5);
plot(pvw(:,1),pvw(:,2),'-k','LineWidth',1.5);
axis equal; axis tight; axis on; colormap jet;

hold on; plot(pvw(:,1),pvw(:,2),'-ko','LineWidth',1.5);

% x = mesh1.dgnodes(:,1,:);
% y = mesh1.dgnodes(:,2,:);
% u = mesh1.dgnodes(:,3,:);
% v = mesh1.dgnodes(:,4,:);
% figure(5); clf; quiver(x,y,u,v); axis equal; axis tight;

scale = 100;
Rout = 45*scale;
n1 = 500;
t1  = linspace(0,2*pi,n1); 
pvc = Rout*[cos(t1(:)) sin(t1(:))];
pvc(:,1) = pvc(:,1) + 5*scale;

ind = find(pvc(:,1)<0);
x0 = pvc(ind,1);
y0 = pvc(ind,2);
figure(1); clf;
plot(x0,y0); 
axis equal; axis tight; axis on; 

ds = 5*ones(2000,1); nref=2;
[x,y,s,Ex,Ey,in,out] = fieldlines(mesh1, mesh1.dgnodes(:,3:4,:), x0, y0, ds, nref, pvc(:,1), pvc(:,2));

ind = 1:1:232;
scale=100;
figure(6); clf; 
hold on;
for i = 1:length(ind)
    int = inpolygon(x(:,ind(i)),y(:,ind(i)),pvc(:,1),pvc(:,2));
    plot(x(int,ind(i))/scale, y(int,ind(i))/scale, '-b', 'LineWidth', 1);
end
fill(pva(:,1)/scale, pva(:,2)/scale,'w');
plot(pva(:,1)/scale, pva(:,2)/scale, '-r', 'LineWidth', 1);
plot(pvw(:,1)/scale, pvw(:,2)/scale, '-r', 'LineWidth', 1);
plot(pvc(:,1)/scale, pvc(:,2)/scale, '-r', 'LineWidth', 1);
hold off;
axis equal; axis tight; axis([-15 25 -20 20]);
set(gca,'FontSize',18);
xlabel('x (cm)','FontSize',20);
ylabel('y (cm)','FontSize',20);


scale = 100;
Rout = 45*scale;
n1 = 160;
t1  = linspace(0.964*pi,0.975*pi,n1); 
pvq = Rout*[cos(t1(:)) sin(t1(:))];
pvq(:,1) = pvq(:,1) + 5*scale;

x1 = pvq(:,1);
y1 = pvq(:,2);
figure(1); clf;
plot(x0,y0,'-r');
hold on;
plot(x1,y1); 
plot(pvw(:,1), pvw(:,2), '-r', 'LineWidth', 5);
plot(pva(:,1), pva(:,2), '-r', 'LineWidth', 1);
axis equal; axis tight; axis on; 

ds = [5*ones(500,1); 2*ones(500,1); 0.5*ones(500,1); 0.25*ones(2000,1); 0.5*ones(500,1); 2*ones(500,1); 5*ones(500,1)]; nref=2;
[x9,y9,s9,Ex,Ey,in,out] = fieldlines(mesh1, mesh1.dgnodes(:,3:4,:), x1, y1, ds, nref, pvc(:,1), pvc(:,2));

scale=100;
ind = 1:1:160;
figure(6); clf; 
hold on;
for i = 1:length(ind)
    %int = inpolygon(x9(:,ind(i)),y9(:,ind(i)),pvc(:,1),pvc(:,2));
    plot(x9(:,ind(i))/scale, y9(:,ind(i))/scale, '-b', 'LineWidth', 1);
end
fill(pva(:,1)/scale, pva(:,2)/scale,'w');
plot(pva(:,1)/scale, pva(:,2)/scale, '-r', 'LineWidth', 1);
plot(pvw(:,1)/scale, pvw(:,2)/scale, '-r', 'LineWidth', 1);
plot(pvc(:,1)/scale, pvc(:,2)/scale, '-r', 'LineWidth', 1);
hold off;
axis equal; axis tight; axis([-15 25 -20 20]);
set(gca,'FontSize',18);
xlabel('x (cm)','FontSize',20);
ylabel('y (cm)','FontSize',20);

return;


