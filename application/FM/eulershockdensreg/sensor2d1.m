function [f,f_udg] = sensor2d1(pg,udg,param,time)
%SENSOR2D1
%    [F,F_UDG] = SENSOR2D1(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    27-Mar-2013 09:10:47
[ng,nc] = size(udg);
nch = 4;
nd = 2;
param1 = param{1};
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
t10 = t5.*t9;
t11 = t2.*t3.*(1.0./2.0);
t12 = t2.*t4.*(1.0./2.0);
t13 = t11+t12;
t24 = t13.*u1;
t14 = -t24+u4;
t15 = param1.*t5.*t6.*t14.*2.0;
t16 = t10+t15;
t17 = param1+1.0;
t18 = 1.0./t17;
t19 = t16.*t18;
t20 = 1.0./sqrt(t19);
t25 = t6.*u2.*u5;
t21 = -t25+u6;
t27 = t6.*u3.*u9;
t22 = -t27+u11;
t23 = 1.0./u1.^3;
t26 = t6.*t21;
t28 = t6.*t22;
t29 = t26+t28;
f = t20.*t29.*1.0e1-5.0;
if nargout > 1
    t30 = 1.0./t19.^(3.0./2.0);
    t31 = t6.*t20.*1.0e1;
    f_udg = [t20.*(t2.*t21+t2.*t22-t23.*u2.*u5-t23.*u3.*u9).*-1.0e1+t18.*t29.*t30.*(t5.*(t3.*t23.*2.0+t4.*t23.*2.0)+param1.*t2.*t5.*t14.*2.0+param1.*t5.*t6.*(t11+t12-u1.*(t3.*t23+t4.*t23)).*2.0).*5.0;t2.*t20.*u5.*-1.0e1-t18.*t29.*t30.*(t2.*t5.*u2.*2.0-param1.*t2.*t5.*u2.*2.0).*5.0;t2.*t20.*u9.*-1.0e1-t18.*t29.*t30.*(t2.*t5.*u3.*2.0-param1.*t2.*t5.*u3.*2.0).*5.0;param1.*t5.*t6.*t18.*t29.*t30.*-1.0e1;t2.*t20.*u2.*-1.0e1;t31;zero;zero;t2.*t20.*u3.*-1.0e1;zero;t31;zero];
end
%f = reshape(f,ng,nch,nd);
%f_udg = reshape(f_udg,ng,nch,nd,nc);