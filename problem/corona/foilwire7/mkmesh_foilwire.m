function [mesh,X1,Y1,X3,Y3] = mkmesh_foilwire(porder,h,xw,yw,rw,scale)

%scale = 1e3/7; % geometry scaling factor
Gamma = 1;
Lambda= 0.1;
delta = 0.1;
k1=1.05;
k2=10;

a  = k1*Gamma*sqrt((1+Lambda)^2 + delta^2); % Inner radius 
b  = k2*Gamma*sqrt((1+Lambda)^2 + delta^2); % Outer radius
c  = Gamma*(-Lambda + 1i*delta);  % origin

r1 = a/a;
r2 = b/a;

n1 = 200;
t1  = linspace(0,2*pi,n1); 
x1 = r1*cos(t1);
y1 = r1*sin(t1);
[X1,Y1] = unitcircle2airfoil(x1,y1,a,real(c),imag(c),scale);

n2 = 90;
t2  = linspace(0,2*pi,n2); 
x2 = r2*cos(t2);
y2 = r2*sin(t2);
[X2,Y2] = unitcircle2airfoil(x2,y2,a,real(c),imag(c),scale);

%xw = 2.0; yw = 0.5; rw = 0.1/4;
n3 = 101;
t3  = linspace(0,2*pi,n3); 
X3 = xw+rw*cos(t3);
Y3 = yw+rw*sin(t3);

pv1 = [X2(:) Y2(:)]; pv1 = pv1(1:end-1,:);
pv2 = [X1(:) Y1(:)]; pv2 = pv2(1:end-1,:);
pv3 = [X3(:) Y3(:)]; pv3 = pv3(1:end-1,:);

[p,t]=polymesh({pv1,pv2,pv3},[1,1,1],[1,0;0,1;0,1],[h,1.2]);
[p,t] = fixmesh(p,t);

s1 = strcat(['all((p(:,1)-' num2str(xw) ').^2+(p(:,2)-' num2str(yw) ').^2<' num2str(rw+1e-3) '^2)']);
s2 = strcat('all((p(:,1)).^2+(p(:,2)).^2>2000^2)');
bndexpr={s1,s2,'true'};   
mesh = mkmesh(p,t,porder,bndexpr,0,1);

% figure(1);clf;
% meshplot(mesh);
% hold on;
% for i = 1:3
%     ind=mesh.f(:,end)==-i;
%     n = mesh.f(ind,1:2);
%     x = mesh.p(unique(n(:)),:);
%     if i==1
%         plot(x(:,1),x(:,2),'or','LineWidth',2);
%     elseif i==2
%         plot(x(:,1),x(:,2),'ob','LineWidth',2);
%     else
%         plot(x(:,1),x(:,2),'og','LineWidth',2);
%     end    
% end


