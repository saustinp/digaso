setapplicationpath('FM/corona');

% See the paper Neimarlija2009.pdf and Medlin1998b.pdf
 
porder = 2;


DV = 40*18.7/12;
Vw = DV+Va;
Rpeek = 4.13*1.2;

% xp=8.5;
% the=30*pi/180;
% R=5.5;
%[mesh1,pva,pvw,xw,yw,rw] = mkmesh_foilwire2(porder,3,xp,the,R);

% stabilization parameter
tau = 4;
eps0 = 1;
nu = 1e-2;
% t1 = 3.4;
% t2 = 5.4;
%t1 = pi/2+pi/4.8;
%t2 = 2*pi+pi/7.0;
% t1 = pi/2+pi/6.0;
% t2 = 2*pi+pi/3.0;
param = {1,nu,eps0,kappa,xw,yw,t1,t2,Wind,tau};

app1.uqpk=0;
hybrid = 'hdg';
app1.source = 'source';
app1.flux = 'flux';
app1.fbou = 'fboucorona';
app1.fhat = 'fhat';
app1.adjoint = 0;
app1.denseblock = 0;
app1.hybrid = hybrid;
app1.localsolve=1;
app1.arg = param;
app1.bcm = [8;1;3];
app1.bcs = [Vw 0;0 0;Va 0];
app1.bcd = [];
app1.bcv = []; 

app1.denseblock = 0;
app1.tdep = false;
app1.wave = false;
app1.alag = false;
app1.flg_q = 1;
app1.flg_p = 0;
app1.flg_g = 0;

app1.fc_q = 1;
app1.fc_u = 0;
app1.fc_p = 0;

app1.np  = 2;
app1.nd  = 2;
app1.nch = 2;                       % Number of componets of UH
app1.nc  = app1.nch*(app1.nd+1);    % Number of componeents of UDG
app1.ncu = 2;

app1.appname = 'corona';
app1.time = [];
app1.dtfc = [];
app1.alpha = [];

% HDG solver
[UDG,UH] = hdg_solve(master,mesh1,app1,UDG,UH,[]);

for i = 1:2
    figure(i); clf; scaplot(mesh1,UDG(:,i,:),[],1); 
    axis equal; axis off; colormap jet; colorbar('FontSize',15);    
    %set(gca, 'LooseInset', get(gca, 'TightInset'));
%     fn = ['numsol0' num2str(i)];
%     print('-dpng',fn);
end

etaw = linspace(0,2*pi,50); 
x0 = cos(etaw)*rw+xw;
y0 = sin(etaw)*rw+yw;
ds = loginc(linspace(0.002,0.1,500),2);
Ua(:,1,:) = UDG(:,1,:);
Ua(:,2,:) = UDG(:,3,:);
Ua(:,3,:) = UDG(:,5,:);    
[Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh1, 0, 1, 0*Ua, Ua, x0, y0, ds, 2, pva(:,1), pva(:,2));        
figure(5);clf; plot(etaw,Peek1,etaw,Rpeek*ones(size(etaw))); 
set(gca,'FontSize',18);
xlabel('\theta','FontSize',20);
ylabel('Peak Integral','FontSize',20);
axis([0 2*pi 3.6 5.1]);

ds = loginc(linspace(0.1,1,700)*5,1);
tetaw = linspace(t1,t2,40); 
x0 = cos(tetaw)*rw+xw;
y0 = sin(tetaw)*rw+yw;
VDG(:,1,:) = UDG(:,3,:)+Wind*mesh1.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,5,:)+Wind*mesh1.dgnodes(:,4,:);
[x,y,s,Ex,Ey,in,out] = fieldlines(mesh1, VDG, x0, y0, ds, 2, pva(:,1), pva(:,2)); 

scale=100;
figure(6); clf; 
plot(x(:,in)/scale, y(:,in)/scale, '-k', 'LineWidth', 1);
hold on;
plot(x(:,out)/scale, y(:,out)/scale, '-b', 'LineWidth', 1);
fill(pva(:,1)/scale, pva(:,2)/scale,'w');
plot(pva(:,1)/scale, pva(:,2)/scale, '-r', 'LineWidth', 1);
plot(pvw(:,1)/scale, pvw(:,2)/scale, '-r', 'LineWidth', 1);
hold off;
axis equal
set(gca,'FontSize',18);
xlabel('x (cm)','FontSize',20);
ylabel('y (cm)','FontSize',20);

return;
