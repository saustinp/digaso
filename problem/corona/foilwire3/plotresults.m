
load windvspotential.mat;
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

figure(1);clf;
plot(Wvel*mu*E0,-Vair*a*E0/1e3,'LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Airfoil Potential (kV)','FontSize',18);
title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;

nref=2;
for i = 1:12
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

curr = zeros(12,1);
for i = 1:12
    fn = ['result' num2str(i) '.mat'];
    load(fn);
    E = sqrt(UDG(:,3,:).^2 + UDG(:,5,:).^2);
    [xdg,u,curr(i)] = surfacedata(mesh,master,E0*mu*UDG(:,2,:).*E,-1);    
end
figure(1);clf;
plot(Wvel*mu*E0,curr/10,'LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Current (\mu A)','FontSize',18);
%title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;

[
    xdg,u] = surfacedata(mesh,master,E,-1);
x = xdg(:,1,:);
y = xdg(:,2,:);
t1 = cart2pol((x-xw),(y-yw));
t1(t1<0) = t1(t1<0)+2*pi;
[t1,ind] = sort(t1(:));
u1 = u(ind); u1 = u1(:);
figure(1);clf;plot(t1(:),u1(:));

