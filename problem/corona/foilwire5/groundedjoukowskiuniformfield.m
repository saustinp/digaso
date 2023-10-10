function [Ex,Ey,u,Ex1,Ey1,u1] = groundedjoukowskiuniformfield(a,Ea,x,y,cx,cy,b,V,scale) 
%a:       Radius of the cylinder (non-dimensional) 
% (cx,cy): center of the cylinder
%Ea:      Uniform electric field pointed in x-direction
%(x, y):  Positions at which the potential and electric field are evaluated
%u:       Potential
%(Ex, Ey):Electric field

% Map Joukowsky airfoil to the unit circle
% [x,y] = inverseJoukowsky(x,y,a,cx,cy);
% x = x-cx;
% y = y-cy;
[x,y,dXdx,dYdx,dXdy,dYdy] = airfoil2unitcircle(x,y,a,cx,cy,scale);
b = b/a;
a = a/a;

r2 = x.^2+y.^2;
% u1 = -Ea*y + (Ea*a^2)*y./r2;
% Ex = (2*Ea*a^2*x.*y)./(r2.^2); 
% Ey = Ea - (Ea*a^2)./r2 + (2*Ea*a^2*y.*y)./(r2.^2); 
u1 = -Ea*x + (Ea*a^2)*x./r2;
Ex = Ea - (Ea*a^2)./r2 + (2*Ea*a^2*x.*x)./(r2.^2);
Ey = (2*Ea*a^2*x.*y)./(r2.^2);  

%Ex1 = Ex; Ey1 = Ey;
Ex1 = Ex.*dXdx + Ey.*dYdx;
Ey1 = Ex.*dXdy + Ey.*dYdy;

if nargin>5
    if (isempty(V)==0) && (isempty(b)==0)
        [ExV,EyV,uV] = groundedcylinderpotential(a,b,V,x,y,0,0); 
        Ex2 = ExV.*dXdx + EyV.*dYdx;
        Ey2 = ExV.*dXdy + EyV.*dYdy;
        Ex = Ex1 + Ex2;
        Ey = Ey1 + Ey2;
        u = uV + u1;
    end
end

Ex = Ex*scale;
Ey = Ey*scale;
u = u*scale;
Ex1 = Ex1*scale;
Ey1 = Ey1*scale;
u1 = u1*scale;


return;



th = 0:0.01:2*pi;
x = cos(th);
y = sin(th);

[Ex,Ey,u] = groundedcylinderuniformfield(1,1,x,y);
figure(1); clf; plot(th, sqrt(Ex.^2+Ey.^2));

En = Ex.*x + Ey.*y;
Et = Ex.*y - Ey.*x;
figure(1); clf; plot(th, Ey);
figure(2); clf; plot(th, Ex);
figure(3); clf; plot(th, En);
figure(4); clf; plot(th, Et);

[Ex,Ey,u] = groundedcylinderuniformfield(1,1,mesh.dgnodes(:,1,:),mesh.dgnodes(:,2,:));

figure(2); clf; scaplot(mesh,Ex,[],1,0);
figure(3); clf; scaplot(mesh,Ey,[],1,0);

syms x y Ea a
r2 = x.^2+y.^2;
u = -Ea*y + (Ea*a^2)*y./r2;
ux = diff(u,'x');
uy = diff(u,'y');
uxx = diff(ux,'x');
uyy = diff(uy,'y');
simplify(uxx+uyy)
