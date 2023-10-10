function [x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref)

% VDG(:,1,:) = UDG(:,3,:) + mesh.dgnodes(:,3,:);
% VDG(:,2,:) = UDG(:,5,:) + mesh.dgnodes(:,4,:);

M = length(x0);
N = length(ds);
x = zeros(M,N+1);
y = zeros(M,N+1);
s = zeros(M,N+1);
Ex = zeros(M,N+1);
Ey = zeros(M,N+1);

x(:,1) = x0(:);
y(:,1) = y0(:);
[udgx] = fieldatx(mesh,VDG,[x(:,1) y(:,1)],nref);
Ex(:,1) = udgx(:,1);
Ey(:,1) = udgx(:,2);

for i = 2:(N+1)       
    Er = sqrt(Ex(:,i-1).^2+Ey(:,i-1).^2);
    x(:,i) = x(:,i-1) + (Ex(:,i-1)./Er)*ds(i-1);
    y(:,i) = y(:,i-1) + (Ey(:,i-1)./Er)*ds(i-1);    
    s(:,i) = s(:,i-1) + sqrt((x(:,i)-x(:,i-1)).^2 + (y(:,i)-y(:,i-1)).^2);
    [udgx] = fieldatx(mesh,VDG,[x(:,i) y(:,i)],nref);
    Ex(:,i) = udgx(:,1);
    Ey(:,i) = udgx(:,2);
end

x = x';
y = y';
s = s';
Ex = Ex';
Ey = Ey';

scale = 1e3/7*2.59; % geometry scaling factor
etaa = linspace(0,2*pi,640); 
Gamma = 1;
Lambda= 0.1;
delta = 0.1;
k1=1.05;
a  = k1*Gamma*sqrt((1+Lambda)^2 + delta^2); % Inner radius 
c  = Gamma*(-Lambda + 1i*delta);  % origin
[xa, ya] = unitcircle2airfoil(cos(etaa),sin(etaa),a,real(c),imag(c),scale);

scale = 1e3/7*2.59; % geometry scaling factor
Rw=550;
xc=-209; yc=60;
xw=xc+Rw*cos(45*pi/180);yw=yc+Rw*sin(45*pi/180); rw = 1;
Nc = 201;
etaw = linspace(0,2*pi,Nc);
xb = cos(etaw)*rw+xw;
yb = sin(etaw)*rw+yw;


in = zeros(M,1);
out =zeros(M,1);
for i = 1:M    
    int = inpolygon(x(:,i),y(:,i),xa,ya);
    if any(int)
        in(i) = 1;
    else
        out(i) = 1;
    end
end
in = find(in==1);
out = find(out==1);

figure(2); clf; 
plot(x(:,in), y(:,in), '-k', 'LineWidth', 1);
hold on;
plot(x(:,out), y(:,out), '-b', 'LineWidth', 1);
fill(xa,ya,'w');
plot(xa, ya, '-r', 'LineWidth', 1);
plot(xb, yb, '-r', 'LineWidth', 1);
hold off;
axis equal


