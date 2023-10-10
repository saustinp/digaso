function [Peek,x,y,s,Ex,Ey] = PeekIntegral(mesh, Vw, Va, UDGw, UDGa, x0, y0, ds, nref, xa, ya)

VDG = Vw*UDGw + Va*UDGa;
[x,y,s,Ex,Ey] = fieldlines(mesh, VDG(:,2:3,:), x0, y0, ds, nref, xa, ya);

e0 = 1;
n = size(s,2);
Peek = zeros(n,1);
for i = 1:n
    e = sqrt(Ex(:,i).^2+Ey(:,i).^2);

    % Find the indices at which e = e0
    ind = bisection(e, e0, 0);

    % Calculate the Peek field using numerical integration
    Peek(i) = numint(e, s(:,i), e0, ind);     
end

return;

R=4.13e-1/0.1; % Threshold for charge emission
scale = 1e3/7*2.59; % geometry scaling factor
Rw=550;
xc=-209; yc=60;
xw=xc+Rw*cos(45*pi/180);yw=yc+Rw*sin(45*pi/180); rw = 1;
N = 101;
etaw = linspace(0,2*pi,N); 
x0 = cos(etaw)*rw+xw;
y0 = sin(etaw)*rw+yw;

%%%%%%%%%%%
DV = 34;
Va = -10;
Vw = DV+Va;
ds = linspace(0.01,0.1,2000)/20;
[Peek,x,y,s,Ex,Ey] = PeekIntegral(mesh, Vw, Va, UDGw, UDGa, x0, y0, ds, 1);
[~,ind] = sort(abs(Peek-R));
figure(4);clf;
plot(etaw,Peek,etaw,R*ones(size(etaw)));
hold on;
plot([etaw(ind(1)) etaw(ind(1))], [min(Peek) R], '-k');
plot([etaw(ind(2)) etaw(ind(2))], [min(Peek) R], '-k');
xlabel('\theta','FontSize',18);
ylabel('Peek Integral','FontSize',18);
title('\Delta V = 34, V_a = -10, V_w = 24','FontSize',18);
set(gca,'FontSize',16);
axis tight;

x=mesh.dgnodes(:,1,:);
y=mesh.dgnodes(:,2,:);
ind = abs(sqrt((x-xw).^2+(y-yw).^2)-rw)<1e-3;
xc=x(ind);
yc=y(ind);
[t1,r1]=cart2pol((xc-xw),(yc-yw));
t1(t1<0) = t1(t1<0)+2*pi;
[t1,ind] = sort(t1);
xc = xc(ind);
yc = yc(ind);
[xc,ind] = unique(xc);
yc = yc(ind);
figure(1);clf;plot(xc,yc,x0,y0);

a1 = linspace(0,2*pi,1000);
b1 = interp1(etaw,Peek,a1);
c1 = b1-R;
c1(c1<0) = 0;
figure(1);clf; plot(etaw,Peek,a1,b1);
save boundata.mat etaw Peek xw yw
figure(1);clf; plot(a1,c1);

ds = linspace(0.1,1,1000)*5;
VDG = Vw*UDGw + Va*UDGa;
[x,y,s,Ex,Ey] = fieldlines(mesh, VDG(:,2:3,:), x0, y0, ds, 1);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
title('\Delta V = 34, V_a = -10, V_w = 24, Wind = 0','FontSize',18);

[Vx, Vy, Psix, Psiy, Psi, Psiax, Psiay, Psia, Vax, Vay] = getvelocity(mesh.dgnodes(:,1,:), mesh.dgnodes(:,2,:), xw, yw, rw, scale);
x=mesh.dgnodes(:,1,:);
[~,imin] = min(x(:));
Vx = Vx/Vx(imin);

Drift(:,1,:) = VDG(:,2,:) + 0.05*Vx;
Drift(:,2,:) = VDG(:,3,:) + 0.05*Vy;
ds = linspace(0.1,1,1000)*2;
[x,y,s,Ex,Ey] = fieldlines(mesh, Drift, x0, y0, ds, 1);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
title('\Delta V = 34, V_a = -10, V_w = 24, Wind = 0.05','FontSize',18);
set(gca,'FontSize',16);

Drift(:,1,:) = VDG(:,2,:) + 0.02*Vx;
Drift(:,2,:) = VDG(:,3,:) + 0.02*Vy;
ds = linspace(0.1,1,1000)*3;
[x,y,s,Ex,Ey] = fieldlines(mesh, Drift, x0, y0, ds, 1);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
title('\Delta V = 34, V_a = -10, V_w = 24, Wind = 0.02','FontSize',18);
set(gca,'FontSize',16);

Drift(:,1,:) = VDG(:,2,:) + 0.01*Vx;
Drift(:,2,:) = VDG(:,3,:) + 0.01*Vy;
ds = linspace(0.1,1,1000)*3;
[x,y,s,Ex,Ey] = fieldlines(mesh, Drift, x0, y0, ds, 1);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
title('\Delta V = 34, V_a = -10, V_w = 24, Wind = 0.01','FontSize',18);
set(gca,'FontSize',16);

Drift(:,1,:) = VDG(:,2,:) + 0.1*Vx;
Drift(:,2,:) = VDG(:,3,:) + 0.1*Vy;
ds = linspace(0.1,1,1000)*2;
[x,y,s,Ex,Ey] = fieldlines(mesh, Drift, x0, y0, ds, 1);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
title('\Delta V = 34, V_a = -10, V_w = 24, Wind = 0.1','FontSize',18);
set(gca,'FontSize',16);


Drift(:,1,:) = VDG(:,2,:) + 0.004*Vx;
Drift(:,2,:) = VDG(:,3,:) + 0.004*Vy;
ds = linspace(0.1,1,1000)*6;
[x,y,s,Ex,Ey] = fieldlines(mesh, Drift, x0, y0, ds, 1);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
title('\Delta V = 34, V_a = -10, V_w = 24, Wind = 0.003','FontSize',18);
set(gca,'FontSize',16);

%%%%%%%%%%%

DV = 34;
Va = -5;
Vw = DV+Va;
ds = linspace(0.01,0.1,1000)/10;
[Peek,x,y,s,Ex,Ey] = PeekIntegral(mesh, Vw, Va, UDGw, UDGa, x0, y0, ds, 1);

figure(1);clf;
plot(etaw,Peek,etaw,R*ones(size(etaw)));
xlabel('\theta','FontSize',18);
ylabel('Peek Integral','FontSize',18);
title('\Delta V = 34, V_a = -5, V_w = 29','FontSize',18);
set(gca,'FontSize',16);
axis tight;

DV = 34;
Va = -15;
Vw = DV+Va;
ds = linspace(0.01,0.1,1000)/10;
[Peek,x,y,s,Ex,Ey] = PeekIntegral(mesh, Vw, Va, UDGw, UDGa, x0, y0, ds, 1);
figure(2);clf;
plot(etaw,Peek,etaw,R*ones(size(etaw)));
xlabel('\theta','FontSize',18);
ylabel('Peek Integral','FontSize',18);
title('\Delta V = 34, V_a = -15, V_w = 19','FontSize',18);
set(gca,'FontSize',16);
axis tight;

UDGp1=U1;
UDGp1(:,2,:)=Ex1;
UDGp1(:,3,:)=Ey1;
UDGp2=U2;
UDGp2(:,2,:)=Ex2;
UDGp2(:,3,:)=Ey2;
[Peek,x,y,s,Ex,Ey] = PeekIntegral(mesh, Vw, Va, UDGw, UDGa, x0, y0, ds, 1);
figure(4);clf;plot(etaw,Peek,etaw,R*ones(size(etaw)));

Ua = Vw*UDGw + Va*UDGa;
Ub = Vw*UDGp1 + Va*UDGp2;

figure(1); clf; scaplot(mesh,Ua(:,1,:),[-10 10],1,0); axis equal; axis tight; colormap jet; axis off;
figure(2); clf; scaplot(mesh,Ub(:,1,:),[-10 10],1,0); axis equal; axis tight; colormap jet; axis off;

figure(1); clf; scaplot(mesh,Ua(:,2,:),[-1 1]/10,1,0); axis equal; axis tight; colormap jet; axis off;
figure(2); clf; scaplot(mesh,Ub(:,2,:),[-1 1]/10,1,0); axis equal; axis tight; colormap jet; axis off;

figure(1); clf; scaplot(mesh,Ua(:,3,:),[-1 1]/10,1,0); axis equal; axis tight; colormap jet; axis off;
figure(2); clf; scaplot(mesh,Ub(:,3,:),[-1 1]/10,1,0); axis equal; axis tight; colormap jet; axis off;


DV = 35;
Vv = -40:5:0;
for i=1:length(Vv)
Va=Vv(i);    
Vw = DV+Va;
ds = linspace(0.01,0.1,2000)/10;
[Peek,x,y,s,Ex,Ey] = PeekIntegral(mesh, Vw, Va, UDGw, UDGa, x0, y0, ds, 1);
[~,ind] = sort(abs(Peek-R));
figure(i);clf;
plot(etaw,Peek,etaw,R*ones(size(etaw)));
% hold on;
% plot([etaw(ind(1)) etaw(ind(1))], [min(Peek) R], '-k');
% plot([etaw(ind(2)) etaw(ind(2))], [min(Peek) R], '-k');
xlabel('\theta','FontSize',18);
ylabel('Peek Integral','FontSize',18);
title(['\Delta V = 35, V_a = ' num2str(Va)],'FontSize',18);
set(gca,'FontSize',16);
axis tight;
end



for i=1:length(Vv)
    figure(i); 
    fn = ['peek1' num2str(i)];
    print('-dpng',fn);
end


Vw = 20;
for i=1:length(Vv)
Va=Vv(i);    
ds = linspace(0.01,0.1,2000)/10;
[Peek,x,y,s,Ex,Ey] = PeekIntegral(mesh, Vw, Va, UDGw, UDGa, x0, y0, ds, 1);
[~,ind] = sort(abs(Peek-R));
figure(i);clf;
plot(etaw,Peek,etaw,R*ones(size(etaw)));
% hold on;
% plot([etaw(ind(1)) etaw(ind(1))], [min(Peek) R], '-k');
% plot([etaw(ind(2)) etaw(ind(2))], [min(Peek) R], '-k');
xlabel('\theta','FontSize',18);
ylabel('Peek Integral','FontSize',18);
title(['V_w = 20, V_a = ' num2str(Va)],'FontSize',18);
set(gca,'FontSize',16);
axis tight;
end

Vw=20;
Wind = [0];



