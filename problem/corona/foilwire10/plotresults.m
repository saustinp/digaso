
load windvspotential.mat;

mu=1.7e-4;
E0=3.0e6;
a=1e-4;
eps0=8.854*1e-12;

figure(1);clf;
plot(Wvel*mu*E0,Vair*a*E0/1e3,'-ob','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Airfoil Potential (-kV)','FontSize',20);
title('\Delta V = 14.8 (kV)','FontSize',20);
set(gca,'FontSize',18);
grid on; box on;
axis([0 50 0 35]);
legend('\theta = 30 deg, R = 5.5 cm','Location','SouthEast');

curr = zeros(13,1);
for i = 1:13
    fn = ['result' num2str(i-1) '.mat'];
    load(fn);
    E = sqrt(UDG(:,3,:).^2 + UDG(:,5,:).^2);
    [xdg,u,curr(i)] = surfacedata(mesh1,master,UDG(:,2,:).*E,-1);    
end
figure(2);clf;
plot(Wvel*mu*E0,curr*eps0*E0*E0*mu*1e6,'-ob','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Current (\mu A)','FontSize',18);
title('\Delta V = 14.8 (kV)','FontSize',18);
set(gca,'FontSize',20);
grid on; box on;
axis([0 50 0 230]);
legend('\theta = 30 deg, R = 5.5 cm','Location','NorthEast');

figure(1);clf;
plot(Wvel1*mu*E0,Vair1*a*E0/1e3,'-ob','LineWidth',1.5);
hold on; plot(Wvel2*mu*E0,Vair2*a*E0/1e3,'-sr','LineWidth',1.5);
hold on; plot(Wvel3*mu*E0,Vair3*a*E0/1e3,'-sr','LineWidth',1.5);
hold on; plot(Wvel4*mu*E0,Vair4*a*E0/1e3,'-sr','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Airfoil Potential (-kV)','FontSize',20);
title('\Delta V = 14.8 (kV)','FontSize',20);
set(gca,'FontSize',18);
grid on; box on;


drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );       

figure(3);clf;
plot(pva(:,1)/scale, pva(:,2)/scale, '-k', 'LineWidth', 1.5);
hold on;
plot(xw/scale, yw/scale, 'ko', 'LineWidth', 1.5, 'MarKerSize', 8);
plot(linspace(0,14,20),0*linspace(0,14,20),'--k','LineWidth', 1.0);
plot(8.5*ones(1,20),linspace(-1,1,20),'--k','LineWidth', 1.0);
% One arrow from left to right with text on left side
xar = [8.5 xw/scale];    % adjust length and location of arrow 
yar = [0.0 yw/scale];      % adjust hieght and width of arrow
drawArrow(xar,yar,'linewidth',2,'color','b');
ang = linspace(0,30*pi/180,100);
xcr = 1.5*cos(ang)+8.5;
ycr = 1.5*sin(ang);
plot(xcr,ycr,'--k','LineWidth', 1.0);
text(10.05,0.5,'\theta','FontSize',18);
text(10.6,1.75,'R','FontSize',18);
text(4,4,'\theta = 30 deg, R = 5.5 cm','FontSize',18);
axis equal
set(gca,'FontSize',18);
xlabel('x (cm)','FontSize',20);
ylabel('y (cm)','FontSize',20);



nref=2;
for i = 1:13
    fn = ['result' num2str(i) '.mat'];
    load(fn);
    ds = loginc(linspace(0.1,2,550)*7,1);
    tetaw = linspace(0,2*pi,200); 
    x0 = cos(tetaw)*rw+xw;
    y0 = sin(tetaw)*rw+yw;
    VDG(:,1,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,3,:);
    VDG(:,2,:) = UDG(:,5,:)+Wind*mesh.dgnodes(:,4,:);
    [x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref);     
    xlabel('x','FontSize',18);
    ylabel('y','FontSize',18);
    title(['Wind = ' num2str(i*0.01)],'FontSize',18);
    set(gca,'FontSize',16);
    axis([-1200 900 -400 1200]);
    fn = ['chargestreamline' num2str(i) '.png'];
    print('-dpng',fn);
end


for i = 1:12
    fn = ['result' num2str(i) '.mat'];
    load(fn);
    etaw = linspace(0,2*pi,100); 
    x0 = cos(etaw)*rw+xw;
    y0 = sin(etaw)*rw+yw;
    nu = 1e-3;
    %ds = loginc(linspace(0.002,0.1,500),2);
    ds = loginc(linspace(0.001,0.12,1000),2);
    Ua(:,1,:) = UDG(:,1,:);
    Ua(:,2,:) = UDG(:,3,:);
    Ua(:,3,:) = UDG(:,5,:);    
    [Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh, 0, 1, 0*Ua, Ua, x0, y0, ds, 2);   
    figure(2);clf;
    plot(etaw,Peek1,etaw,R*ones(size(etaw)));
    xlabel('\theta','FontSize',18);
    ylabel('Peek Integral','FontSize',18);
    title(['Wind = ' num2str(i*0.01)],'FontSize',18);
    set(gca,'FontSize',16);
    axis([0 2*pi 3.85 4.15]);
    fn = ['peekintegral' num2str(i) '.png'];
    print('-dpng',fn);
end


for i = 1:12
    fn = ['result' num2str(i) '.mat'];
    load(fn);
    figure(1); clf; scaplot(mesh,UDG(:,1,:),[Va Vw],2); 
    axis equal; axis([-1000 1200 -800 1200]); axis off; colormap jet; colorbar('FontSize',15);    
    title(['Wind = ' num2str(i*0.01)],'FontSize',18);
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

[xdg,u,uint] = surfacedata(mesh,master,UDG(:,2,:),-1);
x = xdg(:,1,:);
y = xdg(:,2,:);
t1 = cart2pol((x-xw),(y-yw));
t1(t1<0) = t1(t1<0)+2*pi;
[t1,ind] = sort(t1(:));
u1 = u(ind); u1 = u1(:);
figure(1);clf;plot(t1(:),u1(:));
figure(2); clf; scaplot(mesh,UDG(:,2,:),[],1); 
rho = u1;

[xdg,u] = surfacedata(mesh,master,UDG(:,3,:),-1);
x = xdg(:,1,:);
y = xdg(:,2,:);
t1 = cart2pol((x-xw),(y-yw));
t1(t1<0) = t1(t1<0)+2*pi;
[t1,ind] = sort(t1(:));
u1 = u(ind); u1 = u1(:);
figure(1);clf;plot(t1(:),u1(:));
figure(2); clf; scaplot(mesh,UDG(:,3,:),[-0.4 0.4],1); 
ex = u1;

[xdg,u] = surfacedata(mesh,master,UDG(:,5,:),-1);
x = xdg(:,1,:);
y = xdg(:,2,:);
t1 = cart2pol((x-xw),(y-yw));
t1(t1<0) = t1(t1<0)+2*pi;
[t1,ind] = sort(t1(:));
u1 = u(ind); u1 = u1(:);
figure(1);clf;plot(t1(:),u1(:));
figure(2); clf; scaplot(mesh,UDG(:,5,:),[-0.4 0.4],1); 
ey = u1;

x1 = x(:)-xw; x1 = x1(ind);
y1 = y(:)-yw; y1 = y1(ind);
en = ex.*x1 + ey.*y1;
figure(1);clf;plot(t1(:),en(:));

figure(2);clf;plot(t1(:),E0*mu*rho.*en(:));



[
    xdg,u] = surfacedata(mesh,master,E,-1);
x = xdg(:,1,:);
y = xdg(:,2,:);
t1 = cart2pol((x-xw),(y-yw));
t1(t1<0) = t1(t1<0)+2*pi;
[t1,ind] = sort(t1(:));
u1 = u(ind); u1 = u1(:);
figure(1);clf;plot(t1(:),u1(:));

