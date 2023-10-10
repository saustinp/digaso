
load windvspotentialb.mat;
figure(1);clf;
plot(Wvel,Vair);
xlabel('Wind Velocity','FontSize',18);
ylabel('Airfoil Potential','FontSize',18);
title('\Delta V = 40','FontSize',18);
set(gca,'FontSize',16);
axis tight;

mu=1.5e-4;
E0=3.0e6;
a=1e-4;
eps0=8.854*1e-12;


figure(2);clf;
plot(Wvel*mu*E0,-Vair*a*E0/1e3);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Airfoil Potential (kV)','FontSize',18);
title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;

nref=2;
for i = 1:18
    fn = ['result' num2str(i) '.mat'];
    load(fn);
    ds = loginc(linspace(0.1,1,500)*4,1);
    tetaw = linspace(0,2*pi,200); 
    x0 = cos(tetaw)*rw+xw;
    y0 = sin(tetaw)*rw+yw;
    VDG(:,1,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,3,:);
    VDG(:,2,:) = UDG(:,5,:)+Wind*mesh.dgnodes(:,4,:);
    [x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref);     
    xlabel('x','FontSize',18);
    ylabel('y','FontSize',18);
    title(['\Delta V = 40, Wind = ' num2str(i*0.01)],'FontSize',18);
    set(gca,'FontSize',16);
    axis([-800 1050 -600 800]);
    fn = ['chargestreamline' num2str(i) '.png'];
    print('-dpng',fn);
end


for i = 1:18
    fn = ['result' num2str(i) '.mat'];
    load(fn);
    etaw = linspace(0,2*pi,50); 
    x0 = cos(etaw)*rw+xw;
    y0 = sin(etaw)*rw+yw;
    nu = 1e-3;
    ds = loginc(linspace(0.002,0.1,500),2);
    Ua(:,1,:) = UDG(:,1,:);
    Ua(:,2,:) = UDG(:,3,:);
    Ua(:,3,:) = UDG(:,5,:);    
    [Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh, 0, 1, 0*Ua, Ua, x0, y0, ds, 2);   
    figure(2);clf;
    plot(etaw,Peek1,etaw,R*ones(size(etaw)));
    xlabel('\theta','FontSize',18);
    ylabel('Peek Integral','FontSize',18);
    title(['\Delta V = 40, Wind = ' num2str(i*0.01)],'FontSize',18);
    set(gca,'FontSize',16);
    axis([0 2*pi 2.8 4.2]);
    fn = ['peekintegral' num2str(i) '.png'];
    print('-dpng',fn);
end


for i = 1:18
    fn = ['result' num2str(i) '.mat'];
    load(fn);
    figure(1); clf; scaplot(mesh,UDG(:,1,:),[Va Vw],2); 
    axis equal; axis([-1000 1200 -800 1200]); axis off; colormap jet; colorbar('FontSize',15);    
    title(['\Delta V = 40, Wind = ' num2str(i*0.01)],'FontSize',18);
    set(gca,'FontSize',16);
    fn = ['potential' num2str(i) '.png'];
    print('-dpng',fn);
end

for i = 1:18
    fn = ['result' num2str(i) '.mat'];
    load(fn);    
    tm = UDG(:,2,:);
    figure(2); clf; scaplot(mesh,tm,[min(tm(:)) max(tm(:))],2); 
    axis equal; axis([-1000 1200 -800 1200]); axis off; colormap jet; colorbar('FontSize',15);    
    title(['\Delta V = 40, Wind = ' num2str(i*0.01)],'FontSize',18);
    set(gca,'FontSize',16);
    fn = ['chargedensity' num2str(i) '.png'];
    print('-dpng',fn);    
end

curr = zeros(18,1);
for i = 1:18
    fn = ['result' num2str(i) '.mat'];
    load(fn);
    E = sqrt(UDG(:,3,:).^2 + UDG(:,5,:).^2);
    [xdg,u,curr(i)] = surfacedata(mesh,master,E0*mu*UDG(:,2,:).*E,-1);    
end
figure(1);clf;
plot(Wvel*mu*E0,curr,'LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Current (\mu A)','FontSize',18);
%title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;

tetaw = linspace(0,2*pi,200); 
x0 = cos(tetaw)*rw+xw;
y0 = sin(tetaw)*rw+yw;

etaa = linspace(0,2*pi,640); 
Gamma = 1;
Lambda= 0.1;
delta = 0.1;
k1=1.05;
k2=10;
aa  = k1*Gamma*sqrt((1+Lambda)^2 + delta^2); % Inner radius 
b  = k2*Gamma*sqrt((1+Lambda)^2 + delta^2); % Outer radius
c  = Gamma*(-Lambda + 1i*delta);  % origin
[xa, ya] = unitcircle2airfoil(cos(etaa),sin(etaa),aa,real(c),imag(c),scale);

for i = 1:18
fn = ['result' num2str(i) '.mat'];
load(fn);    
tm = UDG(:,2,:)*E0*eps0/a;
figure(2); clf; scaplot(mesh,tm,[0 7e-4],2); 
axis equal; axis([570 810 0 180]); axis off; colormap(flipud(gray)); colorbar('off');
set(gca,'FontSize',16);
hold on;
plot(x0,y0,'-k','LineWidth',1.5);
plot(xa,ya,'-k','LineWidth',1.5);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset*0; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fn = ['chargedensityzoomgray' num2str(i) '.png'];
print('-dpng',fn);    
end


fn = ['result' num2str(i) '.mat'];
load(fn);    
tm = UDG(:,2,:)*E0*eps0/a;
figure(2); clf; scaplot(mesh,tm,[0 7e-4],2); 
axis equal; axis([570 810 0 180]); axis off; colormap(flipud(gray)); colorbar('off');
colorbar('location','NorthOutSide','FontSize',16);  set(gca,'FontSize',16);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset*0; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


