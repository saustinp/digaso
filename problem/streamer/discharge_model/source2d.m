function [s,s_udg] = source2d(p,udg,param,time)
%SOURCE2D
%    [S,S_UDG] = SOURCE2D(P,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    11-Sep-2023 20:24:39
[ng,nc] = size(udg);
nch = 3;
nd = 2;
one = ones(ng,1);
param1 = param{1};
param2 = param{2};
param3 = param{3};
param4 = param{4};
u1 = udg(:,1);
u2 = udg(:,2);
u6 = udg(:,6);
u9 = udg(:,9);
x1 = p(:,1);
zero = zeros(ng,1);
t2 = u6.^2;
t3 = u9.^2;
t4 = 1.0./param1.^2;
t5 = 1.0./param2;
t6 = 1.0./param3;
t9 = param1.*3.4075e+2;
t7 = t6.^3;
t8 = t2+t3;
t11 = one.*param4.*t4.*t6.*x1;
t10 = sqrt(t8);
t12 = 1.0./t10;
t15 = param3.*t10;
t13 = t12.^3;
t14 = t12.^5;
t16 = 1.0./t15.^(1.3e+1./5.0e+1);
t17 = 1.0./t15.^(6.3e+1./5.0e+1);
t18 = t6.*t12.*2.73e+7;
t19 = -t18;
t21 = t7.*t13.*4.3666e+26;
t20 = exp(t19);
t22 = t21+1.1944e+6;
t23 = param1.*t7.*t14.*t20.*u6.*1.30998e+27;
t24 = param1.*t7.*t14.*t20.*u9.*1.30998e+27;
t25 = param1.*t20.*t22;
t26 = -t25;
t28 = t6.*t13.*t25.*u6.*2.73e+7;
t29 = t6.*t13.*t25.*u9.*2.73e+7;
t27 = t9+t26;
t30 = -t28;
t31 = -t29;
t32 = param3.*t5.*t17.*t27.*u1.*u6.*x1.*6.23662e-1;
t33 = param3.*t5.*t17.*t27.*u1.*u9.*x1.*6.23662e-1;
t34 = t5.*t10.*t16.*t27.*x1.*2.3987;
t37 = t5.*t10.*t16.*t27.*u1.*x1.*(-2.3987);
s = [t37;t37;-param4.*t4.*t6.*x1.*(u1-u2)];
if nargout > 1
    t38 = t5.*t12.*t16.*t27.*u1.*u6.*x1.*2.3987;
    t39 = t5.*t12.*t16.*t27.*u1.*u9.*x1.*2.3987;
    t42 = t23+t30;
    t43 = t24+t31;
    t36 = -t34;
    t40 = -t38;
    t41 = -t39;
    t44 = t5.*t10.*t16.*t42.*u1.*x1.*2.3987;
    t45 = t5.*t10.*t16.*t43.*u1.*x1.*2.3987;
    t46 = -t44;
    t47 = -t45;
    t48 = t32+t40+t46;
    t49 = t33+t41+t47;
    s_udg = [t36;t36;-t11;zero;zero;t11;zero;zero;zero;zero;zero;zero;zero;zero;zero;t48;t48;zero;zero;zero;zero;zero;zero;zero;t49;t49;zero];
end
s = reshape(s,ng,nch);
s_udg = reshape(s_udg,ng,nch,nc);