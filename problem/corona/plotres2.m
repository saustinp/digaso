
folders = [5 4 6 13 11 12];
files = [8 7 7 7 5 4];
windall = cell(length(folders),1);
voltall = cell(length(folders),1);
currall = cell(length(folders),1);
for i = 1:length(folders)
    windall{i} = 0:0.01:0.01*10;
    fn = ['foilwire' num2str(folders(i)) '/mesh' num2str(1) '.mat'];
    load(fn);
    for j=1:files(i)+1
        fn = ['foilwire' num2str(folders(i)) '/result' num2str(j-1) '.mat'];
        load(fn);
        voltall{i}(j) = -Va; 
        E = sqrt(UDG(:,3,:).^2 + UDG(:,5,:).^2);
        [~,~,currall{i}(j)] = surfacedata(mesh1,master,UDG(:,2,:).*E,-1);    
    end
    voltall{i}(end:1:11)=voltall{i}(end);
    currall{i}(end:1:11)=0;
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
plot(windall{6}*mu*E0,voltall{6}*a*E0/1e3,'-^k','LineWidth',1.5);
xlabel('Wind Velocity (m/s)','FontSize',20);
ylabel('Airfoil Potential (-kV)','FontSize',20);
set(gca,'FontSize',18);
grid on; box on;
axis([0 45 0 24]);
legend({'\theta = 30, R = 5.5cm', '\theta = 60, R = 5.5cm', '\theta = 90, R = 5.5cm', ...
    '\theta = 10, R = 7cm', '\theta = 30, R = 7cm', '\theta = 60, R = 7cm'},'Location','NorthWest');

the=30*pi/180;
R=5.5;

figure(2);clf;
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



