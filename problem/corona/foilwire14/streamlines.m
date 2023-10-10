function [x,y,s,Ex,Ey] = streamlines(Va, DV, Wind, xw, yw, rw, xq, yq, q, x0, y0, ds, scale)

%[x,y,s,Ex,Ey] = streamlines(-3.8, 68, 0, xw, yw, rw, [], [], [], linspace(1e-1,10,1000), scale);

% Nc = 80;
% etaw = linspace(0,2*pi,Nc);
% x0 = cos(etaw)*rw+xw;
% y0 = sin(etaw)*rw+yw;

Vw = Va+DV;
M = length(x0);
N = length(ds);
x = zeros(M,N+1);
y = zeros(M,N+1);
s = zeros(M,N+1);
Ex = zeros(M,N+1);
Ey = zeros(M,N+1);

x(:,1) = x0(:);
y(:,1) = y0(:);
[Ex(:,1), Ey(:,1)] = getfield(x(:,1), y(:,1), Va, Vw, xw, yw, rw, xq, yq, q, scale);
[Vx(:,1), Vy(:,1)] = getvelocity(x(:,1), y(:,1), xw, yw, rw, scale);
Ex(:,1) = Ex(:,1) + Wind*Vx(:,1);
Ey(:,1) = Ey(:,1) + Wind*Vy(:,1);
for i = 2:(N+1)    
    Er = sqrt(Ex(:,i-1).^2+Ey(:,i-1).^2);
    x(:,i) = x(:,i-1) + (Ex(:,i-1)./Er)*ds(i-1);
    y(:,i) = y(:,i-1) + (Ey(:,i-1)./Er)*ds(i-1);    
    s(:,i) = s(:,i-1) + sqrt((x(:,i)-x(:,i-1)).^2 + (y(:,i)-y(:,i-1)).^2);
    [Ex(:,i), Ey(:,i)] = getfield(x(:,i), y(:,i), Va, Vw, xw, yw, rw, xq, yq, q, scale);
    [Vx(:,i), Vy(:,i)] = getvelocity(x(:,i), y(:,i), xw, yw, rw, scale);
    Ex(:,i) = Ex(:,i) + Wind*Vx(:,i);
    Ey(:,i) = Ey(:,i) + Wind*Vy(:,i);
end

x = x';
y = y';
s = s';
Ex = Ex';
Ey = Ey';

etaa = linspace(0,2*pi,640); 
Gamma = 1;
Lambda= 0.1;
delta = 0.1;
k1=1.05;
k2=10;
a  = k1*Gamma*sqrt((1+Lambda)^2 + delta^2); % Inner radius 
b  = k2*Gamma*sqrt((1+Lambda)^2 + delta^2); % Outer radius
c  = Gamma*(-Lambda + 1i*delta);  % origin
[xa, ya] = unitcircle2airfoil(cos(etaa),sin(etaa),a,real(c),imag(c),scale);

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
hold off;
axis equal



