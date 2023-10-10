

mu=1.5e-4;
E0=3.0e6;
a=1e-4;
eps0=8.854*1e-12;

load data1.mat;

figure(1);clf;
hold on;
plot(windall{1}*mu*E0,voltall{1}*a*E0/1e3,'-ok','LineWidth',1.5);
plot(windall{2}*mu*E0,voltall{2}*a*E0/1e3,'-sk','LineWidth',1.5);
plot(windall{3}*mu*E0,voltall{3}*a*E0/1e3,'->k','LineWidth',1.5);
plot(windall{4}*mu*E0,voltall{4}*a*E0/1e3,'-dk','LineWidth',1.5);
plot(windall{5}*mu*E0,voltall{5}*a*E0/1e3,'-*k','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Airfoil Potential (-kV)','FontSize',20);
set(gca,'FontSize',18);
grid on; box on;
axis([0 85 0 60]);
legend({'13 kV', '14.8 kV', '16.3 kV', '17.5 kV', '18.7 kV'},'Location','NorthWest');

figure(2);clf;
hold on;
plot(windall{1}*mu*E0,currall{1}*eps0*E0*E0*mu*1e6,'-ok','LineWidth',1.5);
plot(windall{2}*mu*E0,currall{2}*eps0*E0*E0*mu*1e6,'-sk','LineWidth',1.5);
plot(windall{3}*mu*E0,currall{3}*eps0*E0*E0*mu*1e6,'->k','LineWidth',1.5);
plot(windall{4}*mu*E0,currall{4}*eps0*E0*E0*mu*1e6,'-dk','LineWidth',1.5);
plot(windall{5}*mu*E0,currall{5}*eps0*E0*E0*mu*1e6,'-*k','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Current (\mu A)','FontSize',18);
set(gca,'FontSize',18);
grid on; box on;
axis([0 85 0 600]);
legend({'13 kV', '14.8 kV', '16.3 kV', '17.5 kV', '18.7 kV'},'Location','NorthEast');

load data2.mat

figure(3);clf;
hold on;
plot(windall{1}*mu*E0,voltall{1}*a*E0/1e3,'-ok','LineWidth',1.5);
plot(windall{2}*mu*E0,voltall{2}*a*E0/1e3,'-sk','LineWidth',1.5);
plot(windall{3}*mu*E0,voltall{3}*a*E0/1e3,'->k','LineWidth',1.5);
plot(windall{4}*mu*E0,voltall{4}*a*E0/1e3,'-dk','LineWidth',1.5);
plot(windall{5}*mu*E0,voltall{5}*a*E0/1e3,'-*k','LineWidth',1.5);
plot(windall{6}*mu*E0,voltall{6}*a*E0/1e3,'-^k','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Airfoil Potential (-kV)','FontSize',20);
set(gca,'FontSize',18);
grid on; box on;
axis([0 45 0 24]);
legend({'\theta = 30, R = 5.5cm', '\theta = 60, R = 5.5cm', '\theta = 90, R = 5.5cm', ...
    '\theta = 10, R = 7cm', '\theta = 30, R = 7cm', '\theta = 60, R = 7cm'},'Location','NorthWest');

figure(4);clf;
hold on;
plot(windall{1}*mu*E0,currall{1}*eps0*E0*E0*mu*1e6,'-ok','LineWidth',1.5);
plot(windall{2}*mu*E0,currall{2}*eps0*E0*E0*mu*1e6,'-sk','LineWidth',1.5);
plot(windall{3}*mu*E0,currall{3}*eps0*E0*E0*mu*1e6,'->k','LineWidth',1.5);
plot(windall{4}*mu*E0,currall{4}*eps0*E0*E0*mu*1e6,'-dk','LineWidth',1.5);
plot(windall{5}*mu*E0,currall{5}*eps0*E0*E0*mu*1e6,'-*k','LineWidth',1.5);
plot(windall{6}*mu*E0,currall{6}*eps0*E0*E0*mu*1e6,'-^k','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Current (\mu A)','FontSize',18);
set(gca,'FontSize',18);
grid on; box on;
axis([0 45 0 110]);
legend({'\theta = 30, R = 5.5cm', '\theta = 60, R = 5.5cm', '\theta = 90, R = 5.5cm', ...
    '\theta = 10, R = 7cm', '\theta = 30, R = 7cm', '\theta = 60, R = 7cm'},'Location','NorthEast');


load data3.mat

figure(5);clf;
hold on;
plot(windall{1}*mu*E0,voltall{1}*a*E0/1e3,'-ok','LineWidth',1.5);
plot(windall{2}*mu*E0,voltall{2}*a*E0/1e3,'-sk','LineWidth',1.5);
plot(windall{3}*mu*E0,voltall{3}*a*E0/1e3,'->k','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Airfoil Potential (-kV)','FontSize',20);
set(gca,'FontSize',18);
grid on; box on;
axis([0 45 0 24]);
legend({'\theta = 150, R = 9cm', '\theta = 170, R = 9cm', '\theta = 180, R = 9cm'},'Location','NorthWest');


figure(6);clf;
hold on;
plot(windall{1}*mu*E0,currall{1}*eps0*E0*E0*mu*1e6,'-ok','LineWidth',1.5);
plot(windall{2}*mu*E0,currall{2}*eps0*E0*E0*mu*1e6,'-sk','LineWidth',1.5);
plot(windall{3}*mu*E0,currall{3}*eps0*E0*E0*mu*1e6,'->k','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Current (\mu A)','FontSize',18);
set(gca,'FontSize',18);
grid on; box on;
axis([0 45 0 110]);
legend({'\theta = 150, R = 9cm', '\theta = 170, R = 9cm', '\theta = 180, R = 9cm'},'Location','NorthEast');



