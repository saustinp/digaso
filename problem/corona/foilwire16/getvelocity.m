function [Vx, Vy, Psix, Psiy, Psi, Psiax, Psiay, Psia, Vax, Vay] = getvelocity(x, y, xw, yw, rw, scale)

alpha = 0;
%scale = 1e3/7; % geometry scaling factor

% Joukowsky profile 
Gamma = 1;
Lambda= 0.1;
delta = 0.1;
k1=1.05;
k2=10;
a  = k1*Gamma*sqrt((1+Lambda)^2 + delta^2); % Inner radius 
b  = k2*Gamma*sqrt((1+Lambda)^2 + delta^2); % Outer radius
c  = Gamma*(-Lambda + 1i*delta);  % origin
cr = real(c); ci = imag(c);

% wire center and radius
%xc = 0.9*scale; yc = 0.3*scale; 
%xw = 2.5*xc; yw = 1.25*yc; rw = 1;
% number of charges on the wire 
nw = 40;

% compute green functions 
[xqw,yqw,xcw,ycw,~,~,uw] = wiregreenfunctionopendomain(a,b,xw,yw,rw,cr,ci,scale,nw,nw) ;

% [xw yw]
% [sqrt((xqw-xw).^2+(yqw-yw).^2)]
% [sqrt((xcw-xw).^2+(ycw-yw).^2)]

% calculate the stream function at the wire due to the airfoil flow 
[~,~,~,ua] = joukowskypotentialflow(xcw,ycw,alpha,a,cr,ci,scale);

% calculate charges at the wire 
%qw = uw\(0.06-ua);
ca = 0.787/1.5; % to make the potential flow around the cylinder 
qw = uw\(ca-ua);

% Electric field due to the airfoil potential
%[Ex,Ey,U] = groundedjoukowskipotential(a,x,y,cr,ci,b,Va,scale);
[Vax,Vay,~,Psi,Psix,Psiy] = joukowskypotentialflow(x,y,alpha,a,cr,ci,scale);
Psia = Psi; Psiax = Psix; Psiay = Psiy;

% Electric field due to the wire
% deltar = 1e-1;
for i = 1:nw    
    [Ewx,Ewy,Uw] = joukowskigreenfunction(a,b,xqw(i),yqw(i),x,y,cr,ci,scale); 
%     ind = ((x-xqw(i)).^2 + (y-yqw(i)).^2)<=deltar^2;    
%     Ewx(ind) = 0;
%     Ewy(ind) = 0;    
    Psi = Psi + Uw*qw(i);
    Psix = Psix + Ewx*qw(i);
    Psiy = Psiy + Ewy*qw(i);
end

Psi = scale*Psi;
Psix = scale*Psix;
Psiy = scale*Psiy;
Psia = scale*Psia;
Psiax = scale*Psiax;
Psiay = scale*Psiay;
Vx = Psiy;
Vy = -Psix;

return;

% position and radius of the wire
porder=1;
scale = 1e3/7*2.59; % geometry scaling factor
Rw=700;
alfa = 160*pi/180;
xc=-209; yc=60;
xw=xc+Rw*cos(alfa);yw=yc+Rw*sin(alfa); rw = 1;
msh = mkmesh_foilwire(porder,400,xw,yw,rw,scale);

x = msh.p(:,1);
y = msh.p(:,2);
[Vx, Vy, Psix, Psiy, Psi, Psiax, Psiay, Psia, Vax, Vay] = getvelocity(x, y, xw, yw, rw, scale);
t = msh.t'; t(4,:) = 1; p = msh.p';
figure(1); clf; pdeplot(p,[],t,'XYData',Psi(:),'Contour','on','Levels',400);  axis tight; axis on; colormap jet;
axis([260 330 -10 50]);

figure(4); clf; pdeplot(p,[],t,'XYData',Psi(:),'FlowData',[Vx(:) Vy(:)],'Contour','on','Levels',400);  axis tight; axis on; colormap jet;

figure(4); clf; pdeplot(p,[],t,'FlowData',[Vx(:) Vy(:)],'Levels',400);  axis tight; axis on; colormap jet;


figure(2); clf; pdeplot(p,[],t,'XYData',Vx(:));  axis tight; axis off; colormap jet; caxis([-2 2]);
axis([260 330 -10 50]);

figure(3); clf; pdeplot(p,[],t,'XYData',Vy(:));  axis tight; axis off; colormap jet; caxis([-2 2]);
axis([260 330 -10 50]);

Vw = -1;      % applied voltage at the wire
Va = -18.5;   % applied voltage at the airfoil
[Ex, Ey, U] = getfield(x, y, Va, Vw, xw, yw, rw, [], [], []);
figure(4); clf; pdeplot(p,[],t,'XYData',U(:),'FlowData',[Ex(:) Ey(:)],'Contour','on','Levels',20);  axis tight; axis on; colormap jet;
figure(5); clf; pdeplot(p,[],t,'XYData',Ex(:));  axis tight; axis off; colormap jet; caxis([-1 1]);
figure(6); clf; pdeplot(p,[],t,'XYData',Ey(:));  axis tight; axis off; colormap jet; caxis([-1 1]);

figure(3); clf; quiver(x,y,Ex,Ey); axis equal; axis tight;

DV = 68;
Va = -10;   % applied voltage at the airfoil
Vw = DV+Va; % applied voltage at the wire
[Ex, Ey, U] = getfield(x, y, Va, Vw, xw, yw, rw, [], [], [], scale);

[Ex1, Ey1, U1] = getfield2(x, y, 1, 0, xw, yw, rw, [], [], [], scale);   
[Ex2, Ey2, U2] = getfield2(x, y, 0, 1, xw, yw, rw, [], [], [], scale);

U3 = Va*U1 + DV*U2;
Ex3 = Va*Ex1 + DV*Ex2;
Ey3 = Va*Ey1 + DV*Ey2;
Wx = Ex3 + 0.0*Vx;
Wy = Ey3 + 0*Vy;
%figure(1);clf;quiver(squeeze(x),squeeze(y),squeeze(Wx),squeeze(Wy));
figure(4); clf; pdeplot(p,[],t,'XYData',U(:),'FlowData',[Wx(:) Wy(:)],'Contour','on','Levels',20);  axis tight; axis on; colormap jet;

figure(1);clf;streamline(squeeze(x),squeeze(y),squeeze(Wx),squeeze(Wy),xw,yw);


x = msh.dgnodes(:,1,:);
y = msh.dgnodes(:,2,:);
[Ex1, Ey1, U1] = getfield2(x, y, 1, 0, xw, yw, rw, [], [], []);   

figure(1); clf; scaplot(msh,U1,[],1,0); axis tight; axis off; colormap jet;
figure(2); clf; scaplot(msh,Ex1,[-0.04 0.04]/4,1,0); axis tight; axis off; colormap jet;
figure(3); clf; scaplot(msh,Ey1,[-0.04 0.04]/4,1,0); axis tight; axis off; colormap jet;


[Ex2, Ey2, U2] = getfield2(x, y, 0, 1, xw, yw, rw, [], [], []);
figure(4); clf; scaplot(msh,U2,[],1,0); axis tight; axis off; colormap jet;
figure(5); clf; scaplot(msh,Ex2,[-0.1 0.1]/5,1,0); axis tight; axis off; colormap jet;
figure(6); clf; scaplot(msh,Ey2,[-0.1 0.1]/5,1,0); axis tight; axis off; colormap jet;


xq = 292; yq = 10; q = 1;
[Ex3, Ey3, U3] = getfield2(x, y, 0, 0, xw, yw, rw, xq, yq, q);
figure(7); clf; scaplot(msh,U3,[0 4],1,0); axis tight; axis off; colormap jet;
figure(8); clf; scaplot(msh,Ex3,[-0.2 0.2],1,0); axis tight; axis off; colormap jet;
figure(9); clf; scaplot(msh,Ey3,[-0.2 0.2],1,0); axis tight; axis off; colormap jet;

load('tm1.mat')
load('tm2.mat')
load('tm3.mat')
ind=100:46000;
ii = I3>0.4; I3(ii)=0.4;
ii = I3<-0.4; I3(ii)=-0.4;
figure(1);clf;plot(time(ind),I3(ind),'LineWidth',1);
axis tight;
xlabel('Time','FontSize',18);
ylabel('Current','FontSize',18);
set(gca,'FontSize',16);


% % position and radius of the wire
% xw = 2; yw = 0.5; rw = 0.025;
% [mesh,xa,ya,xc,yc] = mkmesh_foilwire(1,0.2,xw,yw,rw);
% 
% 
% figure(2); clf; scaplot(mesh,Psiy,[-1 1],1,0); axis tight; axis off; colormap jet;
% figure(3); clf; scaplot(mesh,-Psix,[-1 1],1,0); axis tight; axis off; colormap jet;
% 
% x = mesh.dgnodes(:,1,:);
% y = mesh.dgnodes(:,2,:);
% [Psix, Psiy, Psi, Psiax, Psiay, Psia, Vax, Vay] = getvelocity(x, y, xw, yw, rw);
% master=mkmasterelement(mesh.nd,mesh.porder,mesh.porder,2*mesh.porder,2*mesh.porder,mesh.elemtype,mesh.nodetype);
% 
% 
% qdg = gradu(master.shapnt(:,:,2:end), mesh.dgnodes, Psia);
% figure(1); clf; scaplot(mesh,Psia,[],1,0); axis tight; axis off; colormap jet;
% figure(2); clf; scaplot(mesh,qdg(:,2,:),[-1 1],1,0); axis tight; axis off; colormap jet;
% figure(3); clf; scaplot(mesh,-qdg(:,1,:),[-1 1],1,0); axis tight; axis off; colormap jet;
% figure(4); clf; scaplot(mesh,Vax,[-1 1],1,0); axis tight; axis off; colormap jet;
% figure(5); clf; scaplot(mesh,Vay,[-1 1],1,0); axis tight; axis off; colormap jet;
% 
% qdg = gradu(master.shapnt(:,:,2:end), mesh.dgnodes, Psi);
% figure(1); clf; scaplot(mesh,Psi,[],1,0); axis tight; axis off; colormap jet;
% figure(2); clf; scaplot(mesh,qdg(:,2,:),[-1 1],1,0); axis tight; axis off; colormap jet;
% figure(3); clf; scaplot(mesh,-qdg(:,1,:),[-1 1],1,0); axis tight; axis off; colormap jet;
% figure(4); clf; scaplot(mesh,Psiy,[-1 1],1,0); axis tight; axis off; colormap jet;
% figure(5); clf; scaplot(mesh,-Psix,[-1 1],1,0); axis tight; axis off; colormap jet;
% 
% figure(1); clf; 
% hold on;
% plot(xa,ya,'b-','LineWidth',1.5);
% plot(xc,yc,'b-','LineWidth',1.5);
% quiver(squeeze(x),squeeze(y),squeeze(Psix),squeeze(Psiy));
% %contour(squeeze(x),squeeze(y),squeeze(Psi),150,'LineWidth',1);
% 
% colormap autumn; 
% axis equal;
% axis([-400 800 -400 400]);
% axis off;
% 
% 



