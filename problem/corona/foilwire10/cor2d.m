setapplicationpath('FM/corona');

% See the paper Neimarlija2009.pdf and Medlin1998b.pdf
 
porder = 2;


DV = 40*18.7/12;
Va = -174.65/1.1076;
Vw = DV+Va;
Wind = 0.1;
Rpeek = 4.13*1.2;

xp=8.5;
the=30*pi/180;
R=5.5;
%[mesh1,pva,pvw,xw,yw,rw] = mkmesh_foilwire2(porder,3,xp,the,R);

% stabilization parameter
tau = 4;
eps0 = 1;
kappa = 3.2e-11; %1.25e-2/9;
nu = 1e-2;
% t1 = 3.4;
% t2 = 5.4;
t1 = pi/2+pi/6.0;
t2 = 2*pi+pi/3.0;
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

master = mkmaster(mesh1,2*porder);
[master,mesh1] = preprocess(master,mesh1,hybrid);

UDG = initu(mesh1,{1;0;0;0;0;0});
Ua = Vw*UDGw + Va*UDGa;
UDG(:,1,:) = Ua(:,1,:);
UDG(:,3,:) = Ua(:,2,:);
UDG(:,5,:) = Ua(:,3,:);
UH = inituhat(master,mesh1.elcon,UDG,app1.ncu);

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
axis([0 2*pi 3.6 4.2]);

return;

ds = loginc(linspace(0.1,1,500)*4,1);
tetaw = linspace(t1,t2,100); 
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




% for i = 1:6
%     figure(i); clf; scaplot(mesh1,UDG(:,i,:),[],2); 
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
app1.arg = param;
[UDG,UH] = hdg_solve(master,mesh1,app,UDG,UH,[]);

clm{1} = [];
clm{2} = [];
clm{3} = [-1 1]/5;
clm{5} = [-1 1]/5;
jj = [1 2 3 5];
for i = 1:length(jj)
    figure(i); clf; scaplot(mesh1,UDG(:,jj(i),:),clm{jj(i)},2); 
    axis equal; axis([-800 1600 -800 1200]); axis off; colormap jet; colorbar('FontSize',15);    
    %set(gca, 'LooseInset', get(gca, 'TightInset'));
%     fn = ['numsol0' num2str(i)];
%     print('-dpng',fn);
end
figure(5); clf; scaplot(mesh1,mesh1.dgnodes(:,3,:),[-1 1]/5,2); 
axis equal; axis([-800 1000 -400 800]); axis off; colormap jet; colorbar('FontSize',15);    
figure(6); clf; scaplot(mesh1,mesh1.dgnodes(:,4,:),[-1 1]/5,2); 
axis equal; axis([-800 1000 -400 800]); axis off; colormap jet; colorbar('FontSize',15);    
figure(7); clf; scaplot(mesh1,UDG(:,3,:)+mesh1.dgnodes(:,3,:),[-1 1]/5,2); 
axis equal; axis([-800 1000 -400 800]); axis off; colormap jet; colorbar('FontSize',15);    
figure(8); clf; scaplot(mesh1,UDG(:,5,:)+mesh1.dgnodes(:,4,:),[-1 1]/5,2); 
axis equal; axis([-800 1000 -400 800]); axis off; colormap jet; colorbar('FontSize',15);    


ds = linspace(0.1,1,1000)*2;
VDG(:,1,:) = UDG(:,3,:)+mesh1.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,5,:)+mesh1.dgnodes(:,4,:);
tetaw = linspace(t1,t2,50); 
x1 = cos(tetaw)*rw+xw;
y1 = sin(tetaw)*rw+yw;
[x,y,s,Ex,Ey] = fieldlines(mesh1, VDG, x1, y1, ds, 1);

param = {1,1e-1/100,eps0,kappa,xw,yw,t1,t2,tau};
app1.arg = param;
[UDG,UH] = hdg_solve(master,mesh1,app,UDG,UH,[]);

Ua(:,1,:) = UDG(:,1,:);
Ua(:,2,:) = UDG(:,3,:);
Ua(:,3,:) = UDG(:,5,:);
%ds = linspace(0.01,0.1,5000)/50;
ds = logdec(linspace(1e-3,0.1,200),2);
[Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh1, 0, 1, 0*Ua, Ua, x0, y0, ds, 2);

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
[Peek,x,y,s,Ex,Ey] = PeekIntegral(mesh1, Vw, Va, UDGw, UDGa, x0, y0, ds, 1);
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









