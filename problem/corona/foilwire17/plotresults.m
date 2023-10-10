% volt = zeros(5,1); curr = volt;
% voltb = zeros(5,1); currb = volt;
% files = 5;
% for j=1:files+1
%     fn = ['result' num2str(j-1) '.mat'];
%     load(fn);
%     volt(j) = -Va; 
%     E = sqrt(UDG(:,3,:).^2 + UDG(:,5,:).^2);
%     [~,~,curr(j)] = surfacedata(mesh1,master,UDG(:,2,:).*E,-1);    
%     
%     fn = ['result' num2str(j-1) 'b.mat'];
%     load(fn);
%     voltb(j) = -Va; 
%     E = sqrt(UDG(:,3,:).^2 + UDG(:,5,:).^2);
%     [~,~,currb(j)] = surfacedata(mesh1,master,UDG(:,2,:).*E,-1);    
% end
% % volt(end:1:10)=volt(end);
% % curr(end:1:10)=0;

load R9Theta90.mat;

mu=1.5e-4;
E0=3.0e6;
a=1e-4;
eps0=8.854*1e-12;

Wvel = 0:0.005:0.025;
figure(1);clf;
plot(Wvel*mu*E0,volt*a*E0/1e3,'-ob','LineWidth',1.5); hold on;
plot(Wvel*mu*E0,voltb*a*E0/1e3,'-sk','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Airfoil Potential (-kV)','FontSize',20);
title('\Delta V = 13 (kV)','FontSize',20);
set(gca,'FontSize',18);
grid on; box on;
axis([0 12 -1 6]);
legend('\theta = 90 deg, R = 9 cm (floating)','\theta = 90 deg, R = 9 cm (grounded)', 'Location','NorthWest');


figure(2);clf;
plot(Wvel*mu*E0,curr*eps0*E0*E0*mu*1e6,'-ob','LineWidth',1.5);
hold on; plot(Wvel*mu*E0,currb*eps0*E0*E0*mu*1e6,'-sk','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Current (\mu A)','FontSize',18);
title('\Delta V = 13 (kV)','FontSize',18);
set(gca,'FontSize',20);
grid on; box on;
axis([0 12 0 15]);
legend('\theta = 90 deg, R = 9 cm (floating)','\theta = 90 deg, R = 9 cm (grounded)', 'Location','NorthWest');


