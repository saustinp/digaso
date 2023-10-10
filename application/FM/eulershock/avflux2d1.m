function [f,f_udg] = avflux2d1(pg,udg,param,time)
%AVFLUX2D1
%    [F,F_UDG] = AVFLUX2D1(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    24-Mar-2013 09:51:27
[ng,nc] = size(udg);
nch = 4;
nd = 2;
param1 = param{1};
param8 = param{8};
param9 = param{9};
param10 = param{10};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u4 = udg(:,4);
u5 = udg(:,5);
u6 = udg(:,6);
u9 = udg(:,9);
u11 = udg(:,11);
zero = zeros(ng,1);
t2 = 1.0./u1.^2;
t3 = u2.^2;
t4 = u3.^2;
t5 = param1-1.0;
t6 = 1.0./u1;
t7 = t2.*t3;
t8 = t2.*t4;
t9 = t7+t8;
t10 = t2.*t3.*(1.0./2.0);
t11 = t2.*t4.*(1.0./2.0);
t12 = t10+t11;
t16 = t12.*u1;
t13 = -t16+u4;
t14 = 1.0./param9;
t15 = t5.*t9;
t17 = param1.*t5.*t6.*t13.*2.0;
t18 = t15+t17;
t19 = param1+1.0;
t20 = 1.0./t19;
t21 = t18.*t20;
t22 = 1.0./sqrt(t21);
t34 = t6.*u2.*u5;
t23 = -t34+u6;
t24 = t6.*t23;
t35 = t6.*u3.*u9;
t25 = -t35+u11;
t26 = t6.*t25;
t27 = t24+t26;
t36 = param8.*t22.*t27;
t28 = param10-t36;
t37 = param9.*t28;
t29 = exp(-t37);
t30 = t29+1.0;
t31 = log(t30);
t32 = 1.0./u1.^3;
t33 = param1.*t5.*t6.*t13;
t38 = sqrt(t9);
t39 = sqrt(t33);
t40 = t38+t39;
f = param8.*t14.*t31.*t40;
if nargout > 1
    t41 = t3.*t32.*2.0;
    t42 = t4.*t32.*2.0;
    t43 = t41+t42;
    t44 = t3.*t32;
    t45 = t4.*t32;
    t46 = t44+t45;
    t47 = t10+t11-t46.*u1;
    t48 = 1.0./sqrt(t9);
    t49 = 1.0./sqrt(t33);
    t50 = 1.0./t30;
    t51 = 1.0./t21.^(3.0./2.0);
    t52 = param8.^2;
    t53 = t6.*t22.*t29.*t40.*t50.*t52;
    f_udg = [-param8.*t14.*t31.*(t43.*t48.*(1.0./2.0)+t49.*(param1.*t2.*t5.*t13+param1.*t5.*t6.*t47).*(1.0./2.0))-param8.*t29.*t40.*t50.*(param8.*t22.*(t2.*t23+t2.*t25-t32.*u2.*u5-t32.*u3.*u9)-param8.*t20.*t27.*t51.*(t5.*t43+param1.*t2.*t5.*t13.*2.0+param1.*t5.*t6.*t47.*2.0).*(1.0./2.0));param8.*t14.*t31.*(t2.*t48.*u2-param1.*t2.*t5.*t49.*u2.*(1.0./2.0))-param8.*t29.*t40.*t50.*(param8.*t2.*t22.*u5+param8.*t20.*t27.*t51.*(t2.*t5.*u2.*2.0-param1.*t2.*t5.*u2.*2.0).*(1.0./2.0));param8.*t14.*t31.*(t2.*t48.*u3-param1.*t2.*t5.*t49.*u3.*(1.0./2.0))-param8.*t29.*t40.*t50.*(param8.*t2.*t22.*u9+param8.*t20.*t27.*t51.*(t2.*t5.*u3.*2.0-param1.*t2.*t5.*u3.*2.0).*(1.0./2.0));param1.*param8.*t5.*t6.*t14.*t31.*t49.*(1.0./2.0)-param1.*t5.*t6.*t20.*t27.*t29.*t40.*t50.*t51.*t52;-t2.*t22.*t29.*t40.*t50.*t52.*u2;t53;zero;zero;-t2.*t22.*t29.*t40.*t50.*t52.*u3;zero;t53;zero];
end
f = reshape(f,ng,1);
f_udg = reshape(f_udg,ng,nc);
