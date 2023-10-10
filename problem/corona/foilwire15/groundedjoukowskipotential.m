function [Ex,Ey,u] = groundedjoukowskipotential(a,x,y,cx,cy,b,V,scale) 
%a:       Radius of the cylinder (non-dimensional) 
% (cx,cy): center of the cylinder
%Ea:      Uniform electric field pointed in x-direction
%(x, y):  Positions at which the potential and electric field are evaluated
%u:       Potential
%(Ex, Ey):Electric field

% Map Joukowsky airfoil to the unit circle
[x,y,dXdx,dYdx,dXdy,dYdy] = airfoil2unitcircle(x,y,a,cx,cy,scale);
b = b/a;
a = a/a;

[ExV,EyV,uV] = groundedcylinderpotential(a,b,V,x,y,0,0); 
Ex = ExV.*dXdx + EyV.*dYdx;
Ey = ExV.*dXdy + EyV.*dYdy;
u = uV;


