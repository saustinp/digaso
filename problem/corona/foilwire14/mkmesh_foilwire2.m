function [mesh,pv2,pv3,xw,yw,rw,pv1] = mkmesh_foilwire2(porder,h,xp,the,R)

%x1 = 5.5;
% h = 4;
% xp = 8.5;
% the = 60*pi/180;
% R = 5.5;

scale = 100;
R = scale*R;
xp = scale*xp;
h = scale*h;

Rout = 45*scale;
rw = 1; % 0.1 mm
xw = xp+R*cos(the);
yw = R*sin(the);

% outer circle
n1 = 120;
t1  = linspace(0,2*pi,n1); 
pv1 = Rout*[cos(t1(:)) sin(t1(:))];
pv1(:,1) = pv1(:,1) + 5*scale;
pv1 = pv1(1:end-1,:);
 
% foil
%n2 = 200;
pw2 = [0 0; 0.5 1; 11 1; 14 0; 11 -1; 0.5 -1]*scale; 
t2 = linspace(0,1,10)';
p1 = [pw2(1,1)*(1-t2)+pw2(2,1)*t2  pw2(1,2)*(1-t2) + pw2(2,2)*t2];
t2 = linspace(0,1,80)';
p2 = [pw2(2,1)*(1-t2)+pw2(3,1)*t2  pw2(2,2)*(1-t2) + pw2(3,2)*t2];
t2 = linspace(0,1,30)';
p3 = [pw2(3,1)*(1-t2)+pw2(4,1)*t2  pw2(3,2)*(1-t2) + pw2(4,2)*t2];
t2 = linspace(0,1,30)';
p4 = [pw2(4,1)*(1-t2)+pw2(5,1)*t2  pw2(4,2)*(1-t2) + pw2(5,2)*t2];
t2 = linspace(0,1,60)';
p5 = [pw2(5,1)*(1-t2)+pw2(6,1)*t2  pw2(5,2)*(1-t2) + pw2(6,2)*t2];
t2 = linspace(0,1,10)';
p6 = [pw2(6,1)*(1-t2)+pw2(1,1)*t2  pw2(6,2)*(1-t2) + pw2(1,2)*t2];
pv2 = [p1; p2(2:end,:); p3(2:end,:); p4(2:end,:); p5(2:end,:); p6(2:end-1,:)];
pv2 = pv2(end:-1:1,:);

% wire
n3 = 101;
t3  = linspace(0,2*pi,n3); 
X3 = xw+rw*cos(t3);
Y3 = yw+rw*sin(t3);
pv3 = [X3(:) Y3(:)]; 
pv3 = pv3(1:end-1,:);

figure(1);clf;
plot(pv1(:,1),pv1(:,2),'-b');
hold on;
plot(pv2(:,1),pv2(:,2),'-r');
plot(pv3(:,1),pv3(:,2),'-k');
axis on; axis equal; axis tight;

[p,t]=polymesh({pv1,pv2,pv3},[1,1,1],[1,0;0,1;0,1],[h,1.2]);
[p,t] = fixmesh(p,t);

s1 = strcat(['all((p(:,1)-' num2str(xw) ').^2+(p(:,2)-' num2str(yw) ').^2<' num2str(rw+1e-3) '^2)']);
s2 = strcat('all((p(:,1)).^2+(p(:,2)).^2>2000^2)');
bndexpr={s1,s2,'true'};   
mesh = mkmesh(p,t,porder,bndexpr,0,1);

% figure(2); clf; meshplot(mesh1);

