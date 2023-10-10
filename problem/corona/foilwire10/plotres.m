
clear wind1 wind2 volt1 volt2 curr1 curr2
% ind = 0:1:12;
% for j=1:length(ind)
%     fn = ['result' num2str(ind(j)) '.mat'];
%     load(fn);
%     wind1(j) = Wind;
%     volt1(j) = -Va; 
%     E = sqrt(UDG(:,3,:).^2 + UDG(:,5,:).^2);
%     [~,~,curr1(j)] = surfacedata(mesh1,master,UDG(:,2,:).*E,-1);    
% end
% 
% ind = 0:2:12;
% for j=1:length(ind)
%     fn = ['result' num2str(ind(j)) 'b.mat'];
%     load(fn);
%     wind2(j) = Wind;
%     volt2(j) = -Va; 
%     E = sqrt(UDG(:,3,:).^2 + UDG(:,5,:).^2);
%     [~,~,curr2(j)] = surfacedata(mesh1,master,UDG(:,2,:).*E,-1);    
% end

load data.mat
mu=1.6e-4;
E0=3.0e6;
a=1e-4;
eps0=8.854*1e-12;

figure(1);clf;
hold on;
plot(wind1*mu*E0,volt1*a*E0/1e3,'-sk','LineWidth',1.5,'MarkerSize',8);
plot(wind2*mu*E0,volt2*a*E0/1e3,'-ok','LineWidth',1.5,'MarkerSize',8);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Airfoil Potential (-kV)','FontSize',20);
set(gca,'FontSize',18);
grid on; box on;
axis([0 64 0 50]);
legend({'R = 4.130', 'R = 4.956'},'Location','NorthWest');

figure(2);clf;
hold on;
plot(wind1*mu*E0,curr1*eps0*E0*E0*mu*1e6,'-ok','LineWidth',1.5,'MarkerSize',8);
plot(wind2*mu*E0,curr2*eps0*E0*E0*mu*1e6,'-sk','LineWidth',1.5,'MarkerSize',8);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Current (\mu A)','FontSize',18);
set(gca,'FontSize',18);
grid on; box on;
axis([0 64 0 650]);
legend({'R = 4.130', 'R = 4.956'},'Location','NorthEast');

