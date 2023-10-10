function [Ex,Ey,u] = groundedjoukowskipointcharge(a,b,xc,yc,x,y,cx,cy,scale) 

% Map Joukowsky airfoil to the unit circle
% [x,y] = inverseJoukowsky(x,y,a,cx,cy);
% x = x-cx;
% y = y-cy;
% x = x/a;
% y = y/a;
[x,y,dXdx,dYdx,dXdy,dYdy] = airfoil2unitcircle(x,y,a,cx,cy,scale);
[xc,yc] = airfoil2unitcircle(xc,yc,a,cx,cy,scale);
b = b/a;
a = a/a;


wc = 0;
if wc==0
    [Ex,Ey,u] = groundedcylinderpointcharge1(a,b,xc,yc,x,y);
elseif wc==1
    [Ex,Ey,u] = groundedcylinderpointcharge2(a,b,xc,yc,x,y);    
else
    [Ex,Ey,u] = groundedcylinderpointcharge1(a,b,xc,yc,x,y);
    ind = ((x-xc).^2 + (y-yc).^2)<1e-3;    
    [Ex(ind),Ey(ind),u(ind)] = groundedcylinderpointcharge2(a,b,xc,yc,x(ind),y(ind));    
end

Ex1 = Ex.*dXdx + Ey.*dYdx;
Ey1 = Ex.*dXdy + Ey.*dYdy;
Ex = 2*Ex1;
Ey = 2*Ey1;
u = -u;

% Ex = Ex*scale;
% Ey = Ey*scale;
% u = u*scale;

function [Ex,Ey,u] = groundedcylinderpointcharge1(a,b,xc,yc,x,y)

[t,r] = cart2pol(x,y);
[tc,rc] = cart2pol(xc,yc);

rx = x./r;
ry = y./r;
tx = -y./(r.*r);
ty = x./(r.*r);

h = a/b;
v = exp((t-tc)*1i);
w = r.*v;
v1 = log(w)/(2*1i);
vm = v1 - log(rc)/(2*1i);
vp = v1 + log(rc)/(2*1i);

wx = rx.*v + (w.*tx)*1i;
wy = ry.*v + (w.*ty)*1i;
v3 = (wx./w)*(1/(2*1i));
v4 = (wy./w)*(1/(2*1i));
% vpx = (wx./w)*(1/(2*1i));
% vpy = (wy./w)*(1/(2*1i));

th1m=0; th1mx=0; th1my=0;
th1p=0; th1px=0; th1py=0;
for n=0:2        
    v1 = (2*n+1)*vm;
    v2 = ((-1)^n*2*h^(.25+n*(n+1)))*sin(v1);
    th1m = th1m+v2;
    v2 =((-1)^n*2*h^(.25+n*(n+1))*(2*n+1))*cos(v1);
    th1mx=th1mx + v2.*v3;
    th1my=th1my + v2.*v4;
                
    v1 = (2*n+1)*vp;
    v2 = ((-1)^n*2*h^(.25+n*(n+1)))*sin(v1);
    th1p = th1p + v2;
    v2 = ((-1)^n*2*h^(.25+n*(n+1))*(2*n+1))*cos(v1);
    th1px=th1px + v2.*v3;
    th1py=th1py + v2.*v4;
end

c1 = -(log(rc)/log(h));
v1 = w.^c1;
v2 = th1m./th1p;
v3 = c1*(v1./w);
v4 = th1p.^2;

f = v1.*v2; 
u = 2*log(abs(f));

v1x = wx.*v3;
v1y = wy.*v3;
v2x = (th1mx.*th1p-th1m.*th1px)./v4;
v2y = (th1my.*th1p-th1m.*th1py)./v4;

fx = v1x.*v2 + v1.*v2x;
fy = v1y.*v2 + v1.*v2y;

v1 = real(f).^2 + imag(f).^2;
Ex = (real(fx).*real(f)+imag(fx).*imag(f))./v1;
Ey = (real(fy).*real(f)+imag(fy).*imag(f))./v1;

function [Ex,Ey,u] = groundedcylinderpointcharge2(a,b,xc,yc,x,y)

[t,r] = cart2pol(x,y);
[tc,rc] = cart2pol(xc,yc);

rx = x./r;
ry = y./r;
tx = -y./(r.*r);
ty = x./(r.*r);

rs = min(r,rc);
rl = max(r,rc);

rsx = 0*rs;
rsy = 0*rs;
ind = r<=rc;
rsx(ind) = rx(ind);
rsy(ind) = ry(ind);
rlx = 0*rs;
rly = 0*rs;
ind = r>=rc;
rlx(ind) = rx(ind);
rly(ind) = ry(ind);

a2 = a*a;
b2 = b*b;
u = (1/log(b2/a2))*log(rs.^2./a2).*log(b2./rl.^2);
Ex = (2/log(b2/a2))*((rsx./rs).*log(b2./rl.^2) - log(rs.^2./a2).*(rlx./rl));
Ey = (2/log(b2/a2))*((rsy./rs).*log(b2./rl.^2) - log(rs.^2./a2).*(rly./rl));

M = 50;
% for m = 1:M
%     c = 2/(m*(1 - (a/b)^(2*m)));
%     u1 = cos(m*(t-tc));
%     rsm = rs.^m;
%     u2 = rsm - a^(2*m)./rsm;
%     rlm = rl.^m;
%     u3 = 1./rlm - rlm/(b^(2*m));
%     u = u + c*(u1.*u2.*u3);
%     
%     ua = sin(m*(t-tc));
%     ub = rs.^(m-1) + a^(2*m)./(rs.^(m+1));
%     uc = 1./(rl.^(m+1)) + rl.^(m-1)/(b^(2*m));
%     
%     u1x = -tx.*ua;
%     u2x = rsx.*ub;
%     u3x = -rlx.*uc;
%     Ex = Ex + (m*c)*(u1x.*u2.*u3 + u1.*u2x.*u3 + u1.*u2.*u3x);
%     
%     u1y = -ty.*ua;
%     u2y = rsy.*ub;
%     u3y = -rly.*uc;
%     Ey = Ey + (m*c)*(u1y.*u2.*u3 + u1.*u2y.*u3 + u1.*u2.*u3y);    
% end
rsm = 1;
rlm = 1;
for m = 1:M
    c = 2/(m*(1 - (a/b)^(2*m)));
    
    u0 = m*(t-tc);
    u1 = cos(u0);
    rsm = rsm.*rs;
    u2 = rsm - a^(2*m)./rsm;
    rlm = rlm.*rl;
    u3 = 1./rlm - rlm/(b^(2*m));
    u = u + c*(u1.*u2.*u3);
    
    ua = sin(u0);
    ub = rsm./rs + (a^(2*m))./(rsm.*rs);
    uc = 1./(rlm.*rl) + (rlm./rl)/(b^(2*m));
    
    u1x = -tx.*ua;
    u2x = rsx.*ub;
    u3x = -rlx.*uc;
    u12 = u1.*u2;
    u13 = u1.*u3;
    u23 = u2.*u3;
    Ex = Ex + (m*c)*(u1x.*u23 + u2x.*u13 + u3x.*u12);
    
    u1y = -ty.*ua;
    u2y = rsy.*ub;
    u3y = -rly.*uc;
    Ey = Ey + (m*c)*(u1y.*u23 + u2y.*u13 + u12.*u3y);    
end
Ex = -Ex/2;
Ey = -Ey/2;

return;


rsm = 1;
rlm = 1;
for m = 1:M
    c = 2/(m*(1 - (a/b)^(2*m)));
    
    u0 = m*(t-tc);
    u1 = cos(u0);
    rsm = rsm.*rs;
    u2 = rsm - a^(2*m)./rsm;
    rlm = rlm.*rl;
    u3 = 1./rlm - rlm/(b^(2*m));
    u = u + c*(u1.*u2.*u3);
    
    ua = sin(u0);
    ub = rsm./rs + (a^(2*m))./(rsm.*rs);
    uc = 1./(rlm.*rl) + (rlm./rl)/(b^(2*m));
    
    u1x = -tx.*ua;
    u2x = rsx.*ub;
    u3x = -rlx.*uc;
    u12 = u1.*u2;
    u13 = u1.*u3;
    u23 = u2.*u3;
    Ex = Ex + (m*c)*(u1x.*u23 + u2x.*u13 + u3x.*u12);
    
    u1y = -ty.*ua;
    u2y = rsy.*ub;
    u3y = -rly.*uc;
    Ey = Ey + (m*c)*(u1y.*u23 + u2y.*u13 + u12.*u3y);    
end

a = 1; b = 10; 
rc = 5; thc = pi/4;
xc = rc*cos(thc); yc = rc*sin(thc);

mesh = mkmesh_circincirc(1,101,101,a,b,0.1);
[Ex,Ey,u] = groundedcylinderpointcharge(a,b,xc,yc,mesh.dgnodes(:,1,:),mesh.dgnodes(:,2,:)); 
figure(1); clf; scaplot(mesh,u,[],1,0); colormap(jet);
figure(2); clf; scaplot(mesh,Ex,[-2 2],1,0); colormap(jet);
figure(3); clf; scaplot(mesh,Ey,[-2 2],1,0); colormap(jet);

x=mesh.dgnodes(:,1,:);
y=mesh.dgnodes(:,2,:);
rct = rc;
tct = 0;
dct = tct-thc;
R = [cos(dct) -sin(dct); sin(dct) cos(dct)];
xt = R(1,1)*x+R(1,2)*y;
yt = R(2,1)*x+R(2,2)*y;
mesht = mesh;
mesht.dgnodes(:,1,:) = xt;
mesht.dgnodes(:,2,:) = yt;
figure(4); clf; scaplot(mesht,u,[],1,0); colormap(jet);
figure(5); clf; scaplot(mesht,Ex,[-2 2],1,0); colormap(jet);
figure(6); clf; scaplot(mesht,Ey,[-2 2],1,0); colormap(jet);







master = mkmaster(mesh,2*mesh.porder);
[master,mesh] = preprocess(master,mesh,'hdg');
UHAT=inituhat(master,mesh.elcon,u,1);
QDG = getq(master,mesh,u,UHAT);
figure(4); clf; scaplot(mesh,QDG(:,1,:),[-1 1],1,0);
figure(5); clf; scaplot(mesh,QDG(:,2,:),[-1 1],1,0);


a = 1; b = 100; 
mesh = mkmesh_circincirc(2,101,101,a,b,1);
V = 1;
Ea = 1;
[Ex,Ey,u] = groundedcylinderpotential(a,b,V,mesh.dgnodes(:,1,:),mesh.dgnodes(:,2,:));
Ev = sqrt(Ex.^2+Ey.^2);
[Ex1,Ey1,u] = groundedcylinderuniformfield(a,Ea,mesh.dgnodes(:,1,:),mesh.dgnodes(:,2,:)); 
Ex = Ex1+Ex;
Ey = Ey1+Ey;
figure(1); clf; scaplot(mesh,Ev,[],1,0);
figure(2); clf; scaplot(mesh,Ex,[],1,0);
figure(3); clf; scaplot(mesh,Ey,[],1,0);


x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);

m = 5;
[t,r] = cart2pol(x,y);
r2 = r.^2;
tx = -y./r2;
ty = x./r2;
figure(1); clf; scaplot(mesh,cos(m*t),[],1,0);
figure(2); clf; scaplot(mesh,-m*sin(m*t),[],1,0);
figure(3); clf; scaplot(mesh,m*sin(m*t).*tx,[],1,0);
figure(4); clf; scaplot(mesh,m*sin(m*t).*ty,[],1,0);

master = mkmaster(mesh,2*mesh.porder);
[master,mesh] = preprocess(master,mesh,'hdg');
UHAT=inituhat(master,mesh.elcon,cos(m*t),1);
QDG = getq(master,mesh,cos(m*t),UHAT);
figure(5); clf; scaplot(mesh,QDG(:,1,:),[],1,0);
figure(6); clf; scaplot(mesh,QDG(:,2,:),[],1,0);


a = 1; 
b = 21;
n=1000;
m=1000;
xc = linspace(a,b,n);
x = linspace(a,b,m);
e = x;
for i = i1:i2
    [Ex,Ey] = groundedcylinderpointcharge(a,b,xc(i),0,x,0*x);
    
    Sp=log(xc(i));
    Sm=log(b/xc(i));
    %kp = find(eta>=etashells(n));       
    k = findindex(x, xc(i));
    kp = k:n;
    e(kp) = Sp./(x(kp)*log(b));
    km = 1:(kp(1)-1);            
    e(km) = -Sm./(x(km)*log(b));
    
    figure(1);plot(x,Ex);
    figure(2);plot(x,e);
    pause(0.1);
end


a = 1; 
b = 21;
m=10000;
x = linspace(a,b,m);
e = 0*x; 

xc = 1.01; yc=0;

Sp=log(xc);
Sm=log(b/xc);    
k = findindex(x, xc);
kp = k:m;
e(kp) = Sp./(x(kp)*log(b));
km = 1:(kp(1)-1);            
e(km) = -Sm./(x(km)*log(b));
figure(1);plot(x(1:20),e(1:20));

n = 10;
t = linspace(-pi/10,pi/10,n+1);
t = t(1:end-1);
Ex = 0*x; Ey=0*x;
for i = 1:length(t)
    dct = t(i);
    R = [cos(dct) -sin(dct); sin(dct) cos(dct)];
    xt = R(1,1)*xc+R(1,2)*yc;
    yt = R(2,1)*xc+R(2,2)*yc;
    [Ext,Eyt] = groundedcylinderpointcharge(a,b,xt,yt,x,0*x);    
    Ex = Ex+Ext;
    Ey = Ey+Eyt;
end
figure(2);plot(x(1:20),Ex(1:20)/(n));

figure(3);plot(x,Ey);










