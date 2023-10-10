function [Ex,Ey,u] = groundedcylinderpotential(a,b,V,x,y,cx,cy) 
% a: radius of the cylinder
% (cx,cy): center of the cylinder
% V: Potential at the cylinder
% (x, y):  positions at which the potential and electric field are evaluated
% u: potential
% (Ex, Ey): electric field

% u = A + B*log(r);
% r=a, u = V -> A + B*log(a) = V;
% r=b, u = 0 -> A + B*log(b) = 0;
% B*(log(a)-log(b)) = V -> B = V/log(a/b)
% A = -B*log(b) = -V*log(b)/log(a/b)
% u = B*(log(r) - log(b)) = B*log(r/b) = V*log(r/b)/log(a/b)

r = sqrt((x-cx).^2+(y-cy).^2);
u = (V/log(a/b))*log(r./b);
Ex = -(V/log(a/b))*x./(r.^2);
Ey = -(V/log(a/b))*y./(r.^2);

% Ex = -du/dx = -ds/dy = -(V/log(a/b))*x./((x-cx).^2+(y-cy).^2)
% Ey = -du/dy = ds/dx = -(V/log(a/b))*y./((x-cx).^2+(y-cy).^2)
% s = (V/log(a/b))*x*y./sqrt((x-cx).^2+(y-cy).^2)

