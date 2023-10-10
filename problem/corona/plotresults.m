
load windvspotentialb.mat;
figure(1);clf;
plot(Wvel,Va);
xlabel('Wind Velocity','FontSize',18);
ylabel('Airfoil Potential','FontSize',18);
title('\Delta V = 40','FontSize',18);
set(gca,'FontSize',16);
axis tight;

mu=1.5e-4;
E0=3.0e6;
a=1e-4;

figure(2);clf;
plot(Wvel*mu*E0,Va*a*E0/1e3);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Airfoil Potential (kV)','FontSize',18);
title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;


for i = 1:10
    fn = ['result' num2str(i) 'b.mat'];
    load(fn);
    ds = loginc(linspace(0.1,1,500)*10,1);    
    tetaw = linspace(app.arg{7},app.arg{8},32); 
    x0 = cos(tetaw)*rw+xw;
    y0 = sin(tetaw)*rw+yw;
    VDG(:,1,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,3,:);
    VDG(:,2,:) = UDG(:,5,:)+Wind*mesh.dgnodes(:,4,:);
    [x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref);     
    xlabel('x','FontSize',18);
    ylabel('y','FontSize',18);
    title(['\Delta V = 40, Wind = ' num2str(i*0.01)],'FontSize',18);
    set(gca,'FontSize',16);
    axis([-800 900 -450 850]);
    fn = ['chargestreamline' num2str(i) 'b.png'];
    print('-dpng',fn);
end


for i = 1:10
    fn = ['result' num2str(i) 'b.mat'];
    load(fn);
    etaw = linspace(0,2*pi,50); 
    x0 = cos(etaw)*rw+xw;
    y0 = sin(etaw)*rw+yw;
    ds = linspace(0.01,0.1,2000)/15;
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
    axis([0 2*pi 3.65 4.15]);
    fn = ['peekintegral' num2str(i) 'b.png'];
    print('-dpng',fn);
end


for i = 1:10
    fn = ['result' num2str(i) 'b.mat'];
    load(fn);
    figure(1); clf; scaplot(mesh,UDG(:,1,:),[Va Vw],2); 
    axis equal; axis([-1000 1200 -800 1200]); axis off; colormap jet; colorbar('FontSize',15);    
    title(['\Delta V = 40, Wind = ' num2str(i*0.01)],'FontSize',18);
    set(gca,'FontSize',16);
    fn = ['potential' num2str(i) 'b.png'];
    print('-dpng',fn);
end

for i = 1:10
    fn = ['result' num2str(i) 'b.mat'];
    load(fn);    
    tm = UDG(:,2,:);
    figure(2); clf; scaplot(mesh,tm,[min(tm(:)) max(tm(:))],2); 
    axis equal; axis([-1000 1200 -800 1200]); axis off; colormap jet; colorbar('FontSize',15);    
    title(['\Delta V = 40, Wind = ' num2str(i*0.01)],'FontSize',18);
    set(gca,'FontSize',16);
    fn = ['chargedensity' num2str(i) 'b.png'];
    print('-dpng',fn);    
end

curr = zeros(10,1);
for i = 1:10
    fn = ['result' num2str(i) 'b.mat'];
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


