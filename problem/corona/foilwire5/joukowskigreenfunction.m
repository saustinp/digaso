function [Ex,Ey,u] = joukowskigreenfunction(a,b,xc,yc,x,y,cx,cy,scale) 

[x,y,dXdx,dYdx,dXdy,dYdy] = airfoil2unitcircle(x,y,a,cx,cy,scale);

%figure(1);clf; plot(x(:),y(:),'o'); axis equal; axis tight;

[xc,yc] = airfoil2unitcircle(xc,yc,a,cx,cy,scale);

[Ex,Ey,u] = unitcirclegreenfunction(x,y,xc,yc);

Ex1 = Ex.*dXdx + Ey.*dYdx;
Ey1 = Ex.*dXdy + Ey.*dYdy;
Ex = Ex1;
Ey = Ey1;
%u = -u;

% Ex1 = Gx.*dXdx + Gy.*dYdx;
% Ey1 = Gx.*dXdy + Gy.*dYdy;
% Gx = Ex1;
% Gy = Ey1;



