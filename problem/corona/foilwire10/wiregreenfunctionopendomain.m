function [xc,yc,x,y,ex,ey,uw] = wiregreenfunctionopendomain(a,b,xw,yw,rw,cx,cy,scale,m,n) 

% positions on the circle centered at (xw, yw) and radius rw/2
t = linspace(0,2*pi,m+1);
t = t(1:end-1)';
xc = cos(t)*(rw/2) + xw;
yc = sin(t)*(rw/2) + yw;

% positions on the circle centered at (xw, yw) and radius rw
t = linspace(0,2*pi,n+1);
t = t(1:end-1)';
x = cos(t)*rw + xw;
y = sin(t)*rw + yw;

% potential of the charge at the wire center
ex = zeros(n,m); ey = ex; uw = ex;
for i = 1:m
    [ex(:,i),ey(:,i),uw(:,i)] = joukowskigreenfunction(a,b,xc(i),yc(i),x,y,cx,cy,scale); 
    %[Ex,Ey,u] = joukowskigreenfunction(a,b,xc,yc,x,y,cx,cy,scale) 
end

