
folders = [5 7 8 9 10];
files = [8 10 13 15 16];
windall = cell(length(folders),1);
voltall = cell(length(folders),1);
currall = cell(length(folders),1);
fn = ['foilwire' num2str(folders(1)) '/mesh' num2str(1) '.mat'];
load(fn);
for i = 1:length(folders)
    windall{i} = 0:0.01:0.01*18;
    for j=1:files(i)+1
        fn = ['foilwire' num2str(folders(i)) '/result' num2str(j-1) '.mat'];
        load(fn);
        voltall{i}(j) = -Va; 
        E = sqrt(UDG(:,3,:).^2 + UDG(:,5,:).^2);
        [~,~,currall{i}(j)] = surfacedata(mesh1,master,UDG(:,2,:).*E,-1);    
    end
    voltall{i}(end:1:19)=voltall{i}(end);
    currall{i}(end:1:19)=0;
end

mu=1.5e-4;
E0=3.0e6;
a=1e-4;
eps0=8.854*1e-12;

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


