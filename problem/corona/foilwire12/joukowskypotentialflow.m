function [vx,vy,phi,psi,psix,psiy] = joukowskypotentialflow(x,y,alpha,a,cx,cy,scale)

[x,y,dXdx,dYdx,dXdy,dYdy] = airfoil2unitcircle(x,y,a,cx,cy,scale,alpha);

s = cos(-alpha) + 1i*sin(-alpha);
% z = x + 1i*y;
% g = s*(z + 1./z);
%psi = imag(g);
%phi = real(g);
r2 = x.^2+y.^2;
p = x + x./r2; % potential
q = y - y./r2; % streamline
px = 1 + 1./r2 - (2*x.*x)./(r2.^2);
py = -(2*x.*y)./(r2.^2) ;
qx = (2*x.*y)./(r2.^2) ;
qy = 1 - 1./r2 + (2*y.*y)./(r2.^2);

phi = real(s)*p - imag(s)*q;
phix = real(s)*px - imag(s)*qx;
phiy = real(s)*py - imag(s)*qy;

%vx = phix; vy = phiy;
vx = phix.*dXdx + phiy.*dYdx;
vy = phix.*dXdy + phiy.*dYdy;

psi = real(s)*q - imag(s)*p;
p1x = real(s)*qx - imag(s)*px;
p1y = real(s)*qy - imag(s)*py;

psix = p1x.*dXdx + p1y.*dYdx;
psiy = p1x.*dXdy + p1y.*dYdy;



