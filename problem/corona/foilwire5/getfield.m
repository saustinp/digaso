function [Ex, Ey, U] = getfield(x, y, Va, Vw, xw, yw, rw, xq, yq, q, scale)

%scale = 1e3/7; % geometry scaling factor

% Joukowsky profile 
Gamma = 1;
Lambda= 0.1;
delta = 0.1;
k1=1.05;
k2=10;
a  = k1*Gamma*sqrt((1+Lambda)^2 + delta^2); % Inner radius 
b  = k2*Gamma*sqrt((1+Lambda)^2 + delta^2); % Outer radius
c  = Gamma*(-Lambda + 1i*delta);  % origin
cr = real(c); ci = imag(c);

% wire center and radius
%xc = 0.9*scale; yc = 0.3*scale; 
%xw = 2.5*xc; yw = 1.25*yc; rw = 1;
% number of charges on the wire 
nw = 100;

% compute green functions 
[xqw,yqw,xcw,ycw,~,~,uw] = wiregreenfunction(a,b,xw,yw,rw,cr,ci,scale,nw,nw) ;

% calculate the potential at the wire due to the airfoil applied voltage 
[~,~,ua] = groundedjoukowskipotential(a,xcw,ycw,cr,ci,b,Va,scale);

% calculate the potential at the wire due to the point charges
nq = length(q);
if nq>0
    Uq = zeros(nw,nq);
    for i = 1:nq
        [~,~,Uq(:,i)] = groundedjoukowskipointcharge(a,b,xq(i),yq(i),xcw,ycw,cr,ci,scale); 
    end

    % calculate the charges 
    qw = uw\(Vw-ua-Uq*q(:));
else
    qw = uw\(Vw-ua);
end

% Electric field due to the airfoil potential
[Ex,Ey,U] = groundedjoukowskipotential(a,x,y,cr,ci,b,Va,scale);

% Electric field due to the wire
deltar = 2e-1;
for i = 1:nw
    [Ewx,Ewy,Uw] = groundedjoukowskipointcharge(a,b,xqw(i),yqw(i),x,y,cr,ci,scale); 
    ind = ((x-xqw(i)).^2 + (y-yqw(i)).^2)<=deltar^2;    
    Ewx(ind) = 0;
    Ewy(ind) = 0;    
    U = U + Uw*qw(i);
    Ex = Ex + Ewx*qw(i);
    Ey = Ey + Ewy*qw(i);
end

% Electric field due to the point charges
for i = 1:nq
    [Ewx,Ewy,Uw] = groundedjoukowskipointcharge(a,b,xq(i),yq(i),x,y,cr,ci,scale); 
    ind = ((x-xq(i)).^2 + (y-yq(i)).^2)<=deltar^2;    
    Ewx(ind) = 0;
    Ewy(ind) = 0;    
    U = U + Uw*q(i);
    Ex = Ex + Ewx*q(i);
    Ey = Ey + Ewy*q(i);
end

