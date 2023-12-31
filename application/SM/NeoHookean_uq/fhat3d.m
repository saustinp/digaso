function [fh,fh_udg,fh_uh] = fhat3d(nl,pg,udg,uh,param,time)
%FHAT3D
%    [FH,FH_UDG,FH_UH] = FHAT3D(NL,PG,UDG,UH,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    19-Mar-2017 09:35:51
[ng,nc] = size(udg);
nch = 3;
nd = 3;
lambda = param{2};
mu = param{1};
nl1 = nl(:,1);
nl2 = nl(:,2);
nl3 = nl(:,3);
one = ones(ng,1);
q11 = udg(:,4);
q12 = udg(:,7);
q13 = udg(:,10);
q21 = udg(:,5);
q22 = udg(:,8);
q23 = udg(:,11);
q31 = udg(:,6);
q32 = udg(:,9);
q33 = udg(:,12);
tau = param{3};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
uh1 = uh(:,1);
uh2 = uh(:,2);
uh3 = uh(:,3);
zero = zeros(ng,1);
t2 = q11.*q23.*q32;
t3 = q12.*q21.*q33;
t4 = q13.*q22.*q31;
t6 = q11.*q22.*q33;
t7 = q12.*q23.*q31;
t8 = q13.*q21.*q32;
t5 = t2+t3+t4-t6-t7-t8;
t9 = log(t5);
t12 = lambda.*t9;
t10 = mu-t12;
t11 = 1.0./t5;
t13 = one.*tau;
t14 = q21.*q32;
t18 = q22.*q31;
t15 = t14-t18;
t16 = q22.*q33;
t19 = q23.*q32;
t17 = t16-t19;
t20 = 1.0./(t2+t3+t4-t6-t7-t8).^2;
t21 = q21.*q33;
t23 = q23.*q31;
t22 = t21-t23;
t24 = t17.^2;
t25 = q12.*q33;
t27 = q13.*q32;
t26 = t25-t27;
t28 = q11.*q32;
t30 = q12.*q31;
t29 = t28-t30;
t31 = q11.*q33;
t33 = q13.*q31;
t32 = t31-t33;
t34 = q12.*q23;
t36 = q13.*q22;
t35 = t34-t36;
t37 = q11.*q22;
t39 = q12.*q21;
t38 = t37-t39;
t40 = q11.*q23;
t42 = q13.*q21;
t41 = t40-t42;
fh = [nl1.*(mu.*q11+t10.*t11.*t17)+nl3.*(mu.*q13+t10.*t11.*t15)+nl2.*(mu.*q12-t10.*t11.*t22)+tau.*(u1-uh1);nl1.*(mu.*q21-t10.*t11.*t26)+nl2.*(mu.*q22+t10.*t11.*t32)+nl3.*(mu.*q23-t10.*t11.*t29)+tau.*(u2-uh2);nl1.*(mu.*q31+t10.*t11.*t35)+nl3.*(mu.*q33+t10.*t11.*t38)+nl2.*(mu.*q32-t10.*t11.*t41)+tau.*(u3-uh3)];
if nargout > 1
    t43 = lambda.*t17.*t20.*t26;
    t44 = t10.*t17.*t20.*t26;
    t45 = t43+t44;
    t46 = q32.*t10.*t11;
    t47 = q33.*t10.*t11;
    t48 = t26.^2;
    t49 = lambda.*t17.*t20.*t35;
    t50 = t10.*t17.*t20.*t35;
    t51 = t49+t50;
    t52 = nl1.*t51;
    t53 = q22.*t10.*t11;
    t54 = q23.*t10.*t11;
    t55 = lambda.*t20.*t26.*t35;
    t56 = t10.*t20.*t26.*t35;
    t57 = t55+t56;
    t58 = q12.*t10.*t11;
    t59 = q13.*t10.*t11;
    t60 = t35.^2;
    t61 = lambda.*t17.*t20.*t22;
    t62 = t10.*t17.*t20.*t22;
    t63 = t61+t62;
    t64 = t22.^2;
    t65 = lambda.*t20.*t22.*t26;
    t66 = t10.*t20.*t22.*t26;
    t67 = -t47+t65+t66;
    t68 = lambda.*t20.*t22.*t35;
    t69 = t10.*t20.*t22.*t35;
    t70 = -t54+t68+t69;
    t71 = lambda.*t20.*t22.*t32;
    t72 = t10.*t20.*t22.*t32;
    t73 = t71+t72;
    t74 = q31.*t10.*t11;
    t75 = lambda.*t17.*t20.*t32;
    t76 = t10.*t17.*t20.*t32;
    t77 = t47+t75+t76;
    t78 = lambda.*t20.*t26.*t32;
    t79 = t10.*t20.*t26.*t32;
    t80 = t78+t79;
    t81 = t32.^2;
    t82 = lambda.*t20.*t32.*t35;
    t83 = t10.*t20.*t32.*t35;
    t84 = -t59+t82+t83;
    t85 = lambda.*t20.*t22.*t41;
    t86 = t10.*t20.*t22.*t41;
    t87 = t85+t86;
    t88 = nl2.*t87;
    t89 = q21.*t10.*t11;
    t90 = lambda.*t17.*t20.*t41;
    t91 = t10.*t17.*t20.*t41;
    t92 = t54+t90+t91;
    t93 = lambda.*t20.*t32.*t41;
    t94 = t10.*t20.*t32.*t41;
    t95 = t93+t94;
    t96 = q11.*t10.*t11;
    t97 = lambda.*t20.*t26.*t41;
    t98 = t10.*t20.*t26.*t41;
    t99 = t59+t97+t98;
    t100 = lambda.*t20.*t35.*t41;
    t101 = t10.*t20.*t35.*t41;
    t102 = t100+t101;
    t103 = t41.^2;
    t104 = lambda.*t15.*t20.*t22;
    t105 = t10.*t15.*t20.*t22;
    t106 = t104+t105;
    t107 = lambda.*t15.*t17.*t20;
    t108 = t10.*t15.*t17.*t20;
    t109 = t107+t108;
    t110 = t15.^2;
    t111 = lambda.*t15.*t20.*t32;
    t112 = t10.*t15.*t20.*t32;
    t113 = -t74+t111+t112;
    t114 = lambda.*t15.*t20.*t26;
    t115 = t10.*t15.*t20.*t26;
    t116 = -t46+t114+t115;
    t117 = lambda.*t15.*t20.*t41;
    t118 = t10.*t15.*t20.*t41;
    t119 = -t89+t117+t118;
    t120 = lambda.*t15.*t20.*t35;
    t121 = t10.*t15.*t20.*t35;
    t122 = -t53+t120+t121;
    t123 = lambda.*t15.*t20.*t29;
    t124 = t10.*t15.*t20.*t29;
    t125 = t123+t124;
    t126 = lambda.*t20.*t22.*t29;
    t127 = t10.*t20.*t22.*t29;
    t128 = t74+t126+t127;
    t129 = lambda.*t17.*t20.*t29;
    t130 = t10.*t17.*t20.*t29;
    t131 = t46+t129+t130;
    t132 = lambda.*t20.*t29.*t32;
    t133 = t10.*t20.*t29.*t32;
    t134 = t132+t133;
    t135 = lambda.*t20.*t26.*t29;
    t136 = t10.*t20.*t26.*t29;
    t137 = t135+t136;
    t138 = t29.^2;
    t139 = lambda.*t20.*t29.*t41;
    t140 = t10.*t20.*t29.*t41;
    t141 = -t96+t139+t140;
    t142 = lambda.*t20.*t29.*t35;
    t143 = t10.*t20.*t29.*t35;
    t144 = -t58+t142+t143;
    t145 = lambda.*t15.*t20.*t38;
    t146 = t10.*t15.*t20.*t38;
    t147 = t145+t146;
    t148 = nl3.*t147;
    t149 = lambda.*t20.*t22.*t38;
    t150 = t10.*t20.*t22.*t38;
    t151 = t89+t149+t150;
    t152 = lambda.*t17.*t20.*t38;
    t153 = t10.*t17.*t20.*t38;
    t154 = t53+t152+t153;
    t155 = lambda.*t20.*t29.*t38;
    t156 = t10.*t20.*t29.*t38;
    t157 = t155+t156;
    t158 = lambda.*t20.*t32.*t38;
    t159 = t10.*t20.*t32.*t38;
    t160 = t96+t158+t159;
    t161 = lambda.*t20.*t26.*t38;
    t162 = t10.*t20.*t26.*t38;
    t163 = t58+t161+t162;
    t164 = lambda.*t20.*t38.*t41;
    t165 = t10.*t20.*t38.*t41;
    t166 = t164+t165;
    t167 = lambda.*t20.*t35.*t38;
    t168 = t10.*t20.*t35.*t38;
    t169 = t167+t168;
    t170 = t38.^2;
    fh_udg = [t13;zero;zero;zero;t13;zero;zero;zero;t13;nl1.*(mu+lambda.*t20.*t24+t10.*t20.*t24)-nl2.*t63+nl3.*t109;-nl1.*t45+nl2.*t77-nl3.*t131;t52-nl2.*t92+nl3.*t154;-nl1.*t45+nl2.*t67-nl3.*t116;nl1.*(mu+lambda.*t20.*t48+t10.*t20.*t48)-nl2.*t80+nl3.*t137;-nl1.*t57+nl2.*t99-nl3.*t163;t52-nl2.*t70+nl3.*t122;-nl1.*t57+nl2.*t84-nl3.*t144;nl1.*(mu+lambda.*t20.*t60+t10.*t20.*t60)-nl2.*t102+nl3.*t169;nl2.*(mu+lambda.*t20.*t64+t10.*t20.*t64)-nl1.*t63-nl3.*t106;nl1.*t67-nl2.*t73+nl3.*t128;t88-nl1.*t70-nl3.*t151;-nl2.*t73+nl1.*t77+nl3.*t113;nl2.*(mu+lambda.*t20.*t81+t10.*t20.*t81)-nl1.*t80-nl3.*t134;nl1.*t84-nl2.*t95+nl3.*t160;t88-nl1.*t92-nl3.*t119;-nl2.*t95+nl1.*t99+nl3.*t141;nl2.*(mu+lambda.*t20.*t103+t10.*t20.*t103)-nl1.*t102-nl3.*t166;nl3.*(mu+lambda.*t20.*t110+t10.*t20.*t110)-nl2.*t106+nl1.*t109;nl2.*t113-nl1.*t116-nl3.*t125;t148-nl2.*t119+nl1.*t122;-nl3.*t125+nl2.*t128-nl1.*t131;nl3.*(mu+lambda.*t20.*t138+t10.*t20.*t138)-nl2.*t134+nl1.*t137;nl2.*t141-nl1.*t144-nl3.*t157;t148-nl2.*t151+nl1.*t154;-nl3.*t157+nl2.*t160-nl1.*t163;nl3.*(mu+lambda.*t20.*t170+t10.*t20.*t170)-nl2.*t166+nl1.*t169];
end
if nargout > 2
    fh_uh = [-t13;zero;zero;zero;-t13;zero;zero;zero;-t13];
end
fh = reshape(fh,ng,nch);
fh_udg = reshape(fh_udg,ng,nch,nc);
fh_uh = reshape(fh_uh,ng,nch,nch);
