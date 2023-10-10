
load wireairfoil.mat;
mu=1.5e-4;
E0=3.0e6;
a=1e-4;
eps0=8.854*1e-12;

figure(1);clf;
plot(wind1*mu*E0,Va1*a*E0/1e3);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Airfoil Potential (kV)','FontSize',18);
title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;

figure(2);clf;
plot(wind2*mu*E0,Va2*a*E0/1e3);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Airfoil Potential (kV)','FontSize',18);
title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;

figure(3);clf;
plot(wind3*mu*E0,Va3*a*E0/1e3);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Airfoil Potential (kV)','FontSize',18);
title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;

figure(4);clf;
plot(wind1*mu*E0,current1*eps0*E0*E0*mu*1e6);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Current (\mu A)','FontSize',18);
title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;

figure(5);clf;
plot(wind2*mu*E0,current2*eps0*E0*E0*mu*1e6);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Current (\mu A)','FontSize',18);
title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;

figure(6);clf;
plot(wind3*mu*E0,current3*eps0*E0*E0*mu*1e6);
xlabel('Wind Velocity (m/s)','FontSize',18);
ylabel('Current (\mu A)','FontSize',18);
title('\Delta V = 12 (kV)','FontSize',18);
set(gca,'FontSize',16);
axis tight;


