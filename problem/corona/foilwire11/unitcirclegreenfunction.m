function [Ex,Ey,u] = unitcirclegreenfunction(x,y,xc,yc)

% z = x + 1i*y;
% s = xc + 1i*yc;
% f = (z-s)./(1-z*s);
% G = (1/(2*pi))*log(abs(f));
% 
% fx = - 1./(s*z - 1) - (s*(s - z))./(s*z - 1).^2;
% fy = 1i*fx;
% v1 = real(f).^2 + imag(f).^2;
% Gx = (real(fx).*real(f)+imag(fx).*imag(f))./v1;
% Gy = (real(fy).*real(f)+imag(fy).*imag(f))./v1;
% Gx = (1/(2*pi))*Gx;
% Gy = (1/(2*pi))*Gy;

x0 = xc/(xc^2+yc^2); 
y0 = yc/(xc^2+yc^2);
d1 = sqrt((x-xc).^2+(y-yc).^2);
d2 = sqrt((x-x0).^2+(y-y0).^2);
d3 = sqrt(xc^2+yc^2);
u = (1/(2*pi))*(log(d1)-log(d2)-log(d3));
Ex = (2*x - 2*xc)./(2*((x - xc).^2 + (y - yc).^2)) - (2*x - (2*xc)./(xc.^2 + yc.^2))./(2*((x - xc./(xc.^2 + yc.^2)).^2 + (y - yc./(xc.^2 + yc.^2)).^2));
Ey = (2*y - 2*yc)./(2*((x - xc).^2 + (y - yc).^2)) - (2*y - (2*yc)./(xc.^2 + yc.^2))./(2*((x - xc./(xc.^2 + yc.^2)).^2 + (y - yc./(xc.^2 + yc.^2)).^2));
Ex = (1/(2*pi))*Ex;
Ey = (1/(2*pi))*Ey;

% max(abs(G(:)-u(:)))
% max(abs(Gx(:)-Ex(:)))
% max(abs(Gy(:)-Ey(:)))

