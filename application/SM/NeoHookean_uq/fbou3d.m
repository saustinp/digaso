function [fb,fb_udg,fb_uh] = fbou3d(ib,uinf,nl,pg,udg,uh,param,time)
%FBOU3D.M
%    [FB,FB_UDG,FB_UH] = FBOU3D(IB,UINF,NL,PG,UDG,UH,PARAM,TIME)
[ng,nc] = size(udg);
nch = 3;
nd = 3;
switch (ib)
    case 1
        one = ones(ng,1);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x1 = pg(:,1);
        x2 = pg(:,2);
        x3 = pg(:,3);
        zero = zeros(ng,1);
        fb = [one.*(-uh1+uinf1+x1);-uh2+uinf2+x2;-uh3+uinf3+x3];
        if nargout > 1
            fb_udg = [zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-one;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 2
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
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
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
        fb = [-uinf1+nl1.*(mu.*q11+t10.*t11.*t17)+nl3.*(mu.*q13+t10.*t11.*t15)+nl2.*(mu.*q12-t10.*t11.*t22)+tau.*(u1-uh1);-uinf2+nl1.*(mu.*q21-t10.*t11.*t26)+nl2.*(mu.*q22+t10.*t11.*t32)+nl3.*(mu.*q23-t10.*t11.*t29)+tau.*(u2-uh2);-uinf3+nl1.*(mu.*q31+t10.*t11.*t35)+nl3.*(mu.*q33+t10.*t11.*t38)+nl2.*(mu.*q32-t10.*t11.*t41)+tau.*(u3-uh3)];
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
            fb_udg = [t13;zero;zero;zero;t13;zero;zero;zero;t13;nl1.*(mu+lambda.*t20.*t24+t10.*t20.*t24)-nl2.*t63+nl3.*t109;-nl1.*t45+nl2.*t77-nl3.*t131;t52-nl2.*t92+nl3.*t154;-nl1.*t45+nl2.*t67-nl3.*t116;nl1.*(mu+lambda.*t20.*t48+t10.*t20.*t48)-nl2.*t80+nl3.*t137;-nl1.*t57+nl2.*t99-nl3.*t163;t52-nl2.*t70+nl3.*t122;-nl1.*t57+nl2.*t84-nl3.*t144;nl1.*(mu+lambda.*t20.*t60+t10.*t20.*t60)-nl2.*t102+nl3.*t169;nl2.*(mu+lambda.*t20.*t64+t10.*t20.*t64)-nl1.*t63-nl3.*t106;nl1.*t67-nl2.*t73+nl3.*t128;t88-nl1.*t70-nl3.*t151;-nl2.*t73+nl1.*t77+nl3.*t113;nl2.*(mu+lambda.*t20.*t81+t10.*t20.*t81)-nl1.*t80-nl3.*t134;nl1.*t84-nl2.*t95+nl3.*t160;t88-nl1.*t92-nl3.*t119;-nl2.*t95+nl1.*t99+nl3.*t141;nl2.*(mu+lambda.*t20.*t103+t10.*t20.*t103)-nl1.*t102-nl3.*t166;nl3.*(mu+lambda.*t20.*t110+t10.*t20.*t110)-nl2.*t106+nl1.*t109;nl2.*t113-nl1.*t116-nl3.*t125;t148-nl2.*t119+nl1.*t122;-nl3.*t125+nl2.*t128-nl1.*t131;nl3.*(mu+lambda.*t20.*t138+t10.*t20.*t138)-nl2.*t134+nl1.*t137;nl2.*t141-nl1.*t144-nl3.*t157;t148-nl2.*t151+nl1.*t154;-nl3.*t157+nl2.*t160-nl1.*t163;nl3.*(mu+lambda.*t20.*t170+t10.*t20.*t170)-nl2.*t166+nl1.*t169];
        end
        if nargout > 2
            fb_uh = [-t13;zero;zero;zero;-t13;zero;zero;zero;-t13];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 3
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
        u2 = udg(:,2);
        u3 = udg(:,3);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x1 = pg(:,1);
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
        t14 = q12.*q33;
        t16 = q13.*q32;
        t15 = t14-t16;
        t17 = q22.*q33;
        t22 = q23.*q32;
        t18 = t17-t22;
        t19 = 1.0./(t2+t3+t4-t6-t7-t8).^2;
        t20 = q11.*q32;
        t23 = q12.*q31;
        t21 = t20-t23;
        t24 = q11.*q33;
        t26 = q13.*q31;
        t25 = t24-t26;
        t27 = q12.*q23;
        t29 = q13.*q22;
        t28 = t27-t29;
        t30 = q11.*q22;
        t32 = q12.*q21;
        t31 = t30-t32;
        t33 = q11.*q23;
        t35 = q13.*q21;
        t34 = t33-t35;
        fb = [one.*(-uh1+uinf1+x1);-uinf2+nl1.*(mu.*q21-t10.*t11.*t15)+nl3.*(mu.*q23-t10.*t11.*t21)+nl2.*(mu.*q22+t10.*t11.*t25)+tau.*(u2-uh2);-uinf3+nl1.*(mu.*q31+t10.*t11.*t28)+nl3.*(mu.*q33+t10.*t11.*t31)+nl2.*(mu.*q32-t10.*t11.*t34)+tau.*(u3-uh3)];
        if nargout > 1
            t36 = t15.^2;
            t37 = lambda.*t15.*t19.*t28;
            t38 = t10.*t15.*t19.*t28;
            t39 = t37+t38;
            t40 = q12.*t10.*t11;
            t41 = q13.*t10.*t11;
            t42 = t28.^2;
            t43 = q21.*q33;
            t45 = q23.*q31;
            t44 = t43-t45;
            t46 = q33.*t10.*t11;
            t47 = q23.*t10.*t11;
            t48 = lambda.*t15.*t19.*t25;
            t49 = t10.*t15.*t19.*t25;
            t50 = t48+t49;
            t51 = t25.^2;
            t52 = lambda.*t19.*t25.*t28;
            t53 = t10.*t19.*t25.*t28;
            t54 = -t41+t52+t53;
            t55 = lambda.*t19.*t25.*t34;
            t56 = t10.*t19.*t25.*t34;
            t57 = t55+t56;
            t58 = q11.*t10.*t11;
            t59 = lambda.*t15.*t19.*t34;
            t60 = t10.*t15.*t19.*t34;
            t61 = t41+t59+t60;
            t62 = lambda.*t19.*t28.*t34;
            t63 = t10.*t19.*t28.*t34;
            t64 = t62+t63;
            t65 = t34.^2;
            t66 = q21.*q32;
            t69 = q22.*q31;
            t67 = t66-t69;
            t68 = q31.*t10.*t11;
            t70 = q32.*t10.*t11;
            t71 = q21.*t10.*t11;
            t72 = q22.*t10.*t11;
            t73 = lambda.*t19.*t21.*t25;
            t74 = t10.*t19.*t21.*t25;
            t75 = t73+t74;
            t76 = lambda.*t15.*t19.*t21;
            t77 = t10.*t15.*t19.*t21;
            t78 = t76+t77;
            t79 = t21.^2;
            t80 = lambda.*t19.*t21.*t34;
            t81 = t10.*t19.*t21.*t34;
            t82 = -t58+t80+t81;
            t83 = lambda.*t19.*t21.*t28;
            t84 = t10.*t19.*t21.*t28;
            t85 = -t40+t83+t84;
            t86 = lambda.*t19.*t21.*t31;
            t87 = t10.*t19.*t21.*t31;
            t88 = t86+t87;
            t89 = lambda.*t19.*t25.*t31;
            t90 = t10.*t19.*t25.*t31;
            t91 = t58+t89+t90;
            t92 = lambda.*t15.*t19.*t31;
            t93 = t10.*t15.*t19.*t31;
            t94 = t40+t92+t93;
            t95 = lambda.*t19.*t31.*t34;
            t96 = t10.*t19.*t31.*t34;
            t97 = t95+t96;
            t98 = lambda.*t19.*t28.*t31;
            t99 = t10.*t19.*t28.*t31;
            t100 = t98+t99;
            t101 = t31.^2;
            fb_udg = [zero;zero;zero;zero;t13;zero;zero;zero;t13;zero;-nl1.*(lambda.*t15.*t18.*t19+t10.*t15.*t18.*t19)+nl2.*(t46+lambda.*t18.*t19.*t25+t10.*t18.*t19.*t25)-nl3.*(t70+lambda.*t18.*t19.*t21+t10.*t18.*t19.*t21);nl1.*(lambda.*t18.*t19.*t28+t10.*t18.*t19.*t28)-nl2.*(t47+lambda.*t18.*t19.*t34+t10.*t18.*t19.*t34)+nl3.*(t72+lambda.*t18.*t19.*t31+t10.*t18.*t19.*t31);zero;nl1.*(mu+lambda.*t19.*t36+t10.*t19.*t36)-nl2.*t50+nl3.*t78;-nl1.*t39+nl2.*t61-nl3.*t94;zero;-nl1.*t39+nl2.*t54-nl3.*t85;nl1.*(mu+lambda.*t19.*t42+t10.*t19.*t42)-nl2.*t64+nl3.*t100;zero;nl1.*(-t46+lambda.*t15.*t19.*t44+t10.*t15.*t19.*t44)-nl2.*(lambda.*t19.*t25.*t44+t10.*t19.*t25.*t44)+nl3.*(t68+lambda.*t19.*t21.*t44+t10.*t19.*t21.*t44);-nl1.*(-t47+lambda.*t19.*t28.*t44+t10.*t19.*t28.*t44)+nl2.*(lambda.*t19.*t34.*t44+t10.*t19.*t34.*t44)-nl3.*(t71+lambda.*t19.*t31.*t44+t10.*t19.*t31.*t44);zero;nl2.*(mu+lambda.*t19.*t51+t10.*t19.*t51)-nl1.*t50-nl3.*t75;nl1.*t54-nl2.*t57+nl3.*t91;zero;-nl2.*t57+nl1.*t61+nl3.*t82;nl2.*(mu+lambda.*t19.*t65+t10.*t19.*t65)-nl1.*t64-nl3.*t97;zero;-nl1.*(-t70+lambda.*t15.*t19.*t67+t10.*t15.*t19.*t67)+nl2.*(-t68+lambda.*t19.*t25.*t67+t10.*t19.*t25.*t67)-nl3.*(lambda.*t19.*t21.*t67+t10.*t19.*t21.*t67);nl1.*(-t72+lambda.*t19.*t28.*t67+t10.*t19.*t28.*t67)-nl2.*(-t71+lambda.*t19.*t34.*t67+t10.*t19.*t34.*t67)+nl3.*(lambda.*t19.*t31.*t67+t10.*t19.*t31.*t67);zero;nl3.*(mu+lambda.*t19.*t79+t10.*t19.*t79)-nl2.*t75+nl1.*t78;nl2.*t82-nl1.*t85-nl3.*t88;zero;-nl3.*t88+nl2.*t91-nl1.*t94;nl3.*(mu+lambda.*t19.*t101+t10.*t19.*t101)-nl2.*t97+nl1.*t100];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-t13;zero;zero;zero;-t13];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 4
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
        u3 = udg(:,3);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x2 = pg(:,2);
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
        t25 = q12.*q23;
        t27 = q13.*q22;
        t26 = t25-t27;
        t28 = q11.*q22;
        t30 = q12.*q21;
        t29 = t28-t30;
        t31 = q11.*q23;
        t33 = q13.*q21;
        t32 = t31-t33;
        fb = [-uinf1+nl1.*(mu.*q11+t10.*t11.*t17)+nl3.*(mu.*q13+t10.*t11.*t15)+nl2.*(mu.*q12-t10.*t11.*t22)+tau.*(u1-uh1);-uh2+uinf2+x2;-uinf3+nl1.*(mu.*q31+t10.*t11.*t26)+nl3.*(mu.*q33+t10.*t11.*t29)+nl2.*(mu.*q32-t10.*t11.*t32)+tau.*(u3-uh3)];
        if nargout > 1
            t34 = q12.*q33;
            t36 = q13.*q32;
            t35 = t34-t36;
            t37 = lambda.*t17.*t20.*t26;
            t38 = t10.*t17.*t20.*t26;
            t39 = t37+t38;
            t40 = nl1.*t39;
            t41 = q22.*t10.*t11;
            t42 = q23.*t10.*t11;
            t43 = t26.^2;
            t44 = lambda.*t17.*t20.*t22;
            t45 = t10.*t17.*t20.*t22;
            t46 = t44+t45;
            t47 = t22.^2;
            t48 = lambda.*t20.*t22.*t26;
            t49 = t10.*t20.*t22.*t26;
            t50 = -t42+t48+t49;
            t51 = q11.*q33;
            t53 = q13.*q31;
            t52 = t51-t53;
            t54 = q13.*t10.*t11;
            t55 = lambda.*t20.*t22.*t32;
            t56 = t10.*t20.*t22.*t32;
            t57 = t55+t56;
            t58 = nl2.*t57;
            t59 = q21.*t10.*t11;
            t60 = lambda.*t17.*t20.*t32;
            t61 = t10.*t17.*t20.*t32;
            t62 = t42+t60+t61;
            t63 = lambda.*t20.*t26.*t32;
            t64 = t10.*t20.*t26.*t32;
            t65 = t63+t64;
            t66 = t32.^2;
            t67 = lambda.*t15.*t20.*t22;
            t68 = t10.*t15.*t20.*t22;
            t69 = t67+t68;
            t70 = lambda.*t15.*t17.*t20;
            t71 = t10.*t15.*t17.*t20;
            t72 = t70+t71;
            t73 = t15.^2;
            t74 = lambda.*t15.*t20.*t32;
            t75 = t10.*t15.*t20.*t32;
            t76 = -t59+t74+t75;
            t77 = lambda.*t15.*t20.*t26;
            t78 = t10.*t15.*t20.*t26;
            t79 = -t41+t77+t78;
            t80 = q11.*q32;
            t82 = q12.*q31;
            t81 = t80-t82;
            t83 = q11.*t10.*t11;
            t84 = q12.*t10.*t11;
            t85 = lambda.*t15.*t20.*t29;
            t86 = t10.*t15.*t20.*t29;
            t87 = t85+t86;
            t88 = nl3.*t87;
            t89 = lambda.*t20.*t22.*t29;
            t90 = t10.*t20.*t22.*t29;
            t91 = t59+t89+t90;
            t92 = lambda.*t17.*t20.*t29;
            t93 = t10.*t17.*t20.*t29;
            t94 = t41+t92+t93;
            t95 = lambda.*t20.*t29.*t32;
            t96 = t10.*t20.*t29.*t32;
            t97 = t95+t96;
            t98 = lambda.*t20.*t26.*t29;
            t99 = t10.*t20.*t26.*t29;
            t100 = t98+t99;
            t101 = t29.^2;
            fb_udg = [t13;zero;zero;zero;zero;zero;zero;zero;t13;nl1.*(mu+lambda.*t20.*t24+t10.*t20.*t24)-nl2.*t46+nl3.*t72;zero;t40-nl2.*t62+nl3.*t94;-nl3.*(-q32.*t10.*t11+lambda.*t15.*t20.*t35+t10.*t15.*t20.*t35)+nl2.*(-q33.*t10.*t11+lambda.*t20.*t22.*t35+t10.*t20.*t22.*t35)-nl1.*(lambda.*t17.*t20.*t35+t10.*t17.*t20.*t35);zero;-nl1.*(lambda.*t20.*t26.*t35+t10.*t20.*t26.*t35)+nl2.*(t54+lambda.*t20.*t32.*t35+t10.*t20.*t32.*t35)-nl3.*(t84+lambda.*t20.*t29.*t35+t10.*t20.*t29.*t35);t40-nl2.*t50+nl3.*t79;zero;nl1.*(mu+lambda.*t20.*t43+t10.*t20.*t43)-nl2.*t65+nl3.*t100;nl2.*(mu+lambda.*t20.*t47+t10.*t20.*t47)-nl1.*t46-nl3.*t69;zero;t58-nl1.*t50-nl3.*t91;nl3.*(-q31.*t10.*t11+lambda.*t15.*t20.*t52+t10.*t15.*t20.*t52)+nl1.*(q33.*t10.*t11+lambda.*t17.*t20.*t52+t10.*t17.*t20.*t52)-nl2.*(lambda.*t20.*t22.*t52+t10.*t20.*t22.*t52);zero;nl1.*(-t54+lambda.*t20.*t26.*t52+t10.*t20.*t26.*t52)-nl2.*(lambda.*t20.*t32.*t52+t10.*t20.*t32.*t52)+nl3.*(t83+lambda.*t20.*t29.*t52+t10.*t20.*t29.*t52);t58-nl1.*t62-nl3.*t76;zero;nl2.*(mu+lambda.*t20.*t66+t10.*t20.*t66)-nl1.*t65-nl3.*t97;nl3.*(mu+lambda.*t20.*t73+t10.*t20.*t73)-nl2.*t69+nl1.*t72;zero;t88-nl2.*t76+nl1.*t79;-nl1.*(q32.*t10.*t11+lambda.*t17.*t20.*t81+t10.*t17.*t20.*t81)+nl2.*(q31.*t10.*t11+lambda.*t20.*t22.*t81+t10.*t20.*t22.*t81)-nl3.*(lambda.*t15.*t20.*t81+t10.*t15.*t20.*t81);zero;-nl1.*(-t84+lambda.*t20.*t26.*t81+t10.*t20.*t26.*t81)+nl2.*(-t83+lambda.*t20.*t32.*t81+t10.*t20.*t32.*t81)-nl3.*(lambda.*t20.*t29.*t81+t10.*t20.*t29.*t81);t88-nl2.*t91+nl1.*t94;zero;nl3.*(mu+lambda.*t20.*t101+t10.*t20.*t101)-nl2.*t97+nl1.*t100];
        end
        if nargout > 2
            fb_uh = [-t13;zero;zero;zero;-one;zero;zero;zero;-t13];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 5
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
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x3 = pg(:,3);
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
        fb = [-uinf1+nl1.*(mu.*q11+t10.*t11.*t17)+nl3.*(mu.*q13+t10.*t11.*t15)+nl2.*(mu.*q12-t10.*t11.*t22)+tau.*(u1-uh1);-uinf2+nl1.*(mu.*q21-t10.*t11.*t26)+nl2.*(mu.*q22+t10.*t11.*t32)+nl3.*(mu.*q23-t10.*t11.*t29)+tau.*(u2-uh2);-uh3+uinf3+x3];
        if nargout > 1
            t34 = lambda.*t17.*t20.*t26;
            t35 = t10.*t17.*t20.*t26;
            t36 = t34+t35;
            t37 = q32.*t10.*t11;
            t38 = q33.*t10.*t11;
            t39 = t26.^2;
            t40 = q12.*q23;
            t42 = q13.*q22;
            t41 = t40-t42;
            t43 = lambda.*t17.*t20.*t22;
            t44 = t10.*t17.*t20.*t22;
            t45 = t43+t44;
            t46 = t22.^2;
            t47 = lambda.*t20.*t22.*t26;
            t48 = t10.*t20.*t22.*t26;
            t49 = -t38+t47+t48;
            t50 = lambda.*t20.*t22.*t32;
            t51 = t10.*t20.*t22.*t32;
            t52 = t50+t51;
            t53 = q31.*t10.*t11;
            t54 = lambda.*t17.*t20.*t32;
            t55 = t10.*t17.*t20.*t32;
            t56 = t38+t54+t55;
            t57 = lambda.*t20.*t26.*t32;
            t58 = t10.*t20.*t26.*t32;
            t59 = t57+t58;
            t60 = t32.^2;
            t61 = q11.*q23;
            t63 = q13.*q21;
            t62 = t61-t63;
            t64 = lambda.*t15.*t20.*t22;
            t65 = t10.*t15.*t20.*t22;
            t66 = t64+t65;
            t67 = lambda.*t15.*t17.*t20;
            t68 = t10.*t15.*t17.*t20;
            t69 = t67+t68;
            t70 = t15.^2;
            t71 = lambda.*t15.*t20.*t32;
            t72 = t10.*t15.*t20.*t32;
            t73 = -t53+t71+t72;
            t74 = lambda.*t15.*t20.*t26;
            t75 = t10.*t15.*t20.*t26;
            t76 = -t37+t74+t75;
            t77 = lambda.*t15.*t20.*t29;
            t78 = t10.*t15.*t20.*t29;
            t79 = t77+t78;
            t80 = lambda.*t20.*t22.*t29;
            t81 = t10.*t20.*t22.*t29;
            t82 = t53+t80+t81;
            t83 = lambda.*t17.*t20.*t29;
            t84 = t10.*t17.*t20.*t29;
            t85 = t37+t83+t84;
            t86 = lambda.*t20.*t29.*t32;
            t87 = t10.*t20.*t29.*t32;
            t88 = t86+t87;
            t89 = lambda.*t20.*t26.*t29;
            t90 = t10.*t20.*t26.*t29;
            t91 = t89+t90;
            t92 = t29.^2;
            t93 = q11.*q22;
            t95 = q12.*q21;
            t94 = t93-t95;
            fb_udg = [t13;zero;zero;zero;t13;zero;zero;zero;zero;nl1.*(mu+lambda.*t20.*t24+t10.*t20.*t24)-nl2.*t45+nl3.*t69;-nl1.*t36+nl2.*t56-nl3.*t85;zero;-nl1.*t36+nl2.*t49-nl3.*t76;nl1.*(mu+lambda.*t20.*t39+t10.*t20.*t39)-nl2.*t59+nl3.*t91;zero;nl3.*(-q22.*t10.*t11+lambda.*t15.*t20.*t41+t10.*t15.*t20.*t41)-nl2.*(-q23.*t10.*t11+lambda.*t20.*t22.*t41+t10.*t20.*t22.*t41)+nl1.*(lambda.*t17.*t20.*t41+t10.*t17.*t20.*t41);-nl3.*(-q12.*t10.*t11+lambda.*t20.*t29.*t41+t10.*t20.*t29.*t41)+nl2.*(-q13.*t10.*t11+lambda.*t20.*t32.*t41+t10.*t20.*t32.*t41)-nl1.*(lambda.*t20.*t26.*t41+t10.*t20.*t26.*t41);zero;nl2.*(mu+lambda.*t20.*t46+t10.*t20.*t46)-nl1.*t45-nl3.*t66;nl1.*t49-nl2.*t52+nl3.*t82;zero;-nl2.*t52+nl1.*t56+nl3.*t73;nl2.*(mu+lambda.*t20.*t60+t10.*t20.*t60)-nl1.*t59-nl3.*t88;zero;-nl3.*(-q21.*t10.*t11+lambda.*t15.*t20.*t62+t10.*t15.*t20.*t62)-nl1.*(q23.*t10.*t11+lambda.*t17.*t20.*t62+t10.*t17.*t20.*t62)+nl2.*(lambda.*t20.*t22.*t62+t10.*t20.*t22.*t62);nl1.*(q13.*t10.*t11+lambda.*t20.*t26.*t62+t10.*t20.*t26.*t62)+nl3.*(-q11.*t10.*t11+lambda.*t20.*t29.*t62+t10.*t20.*t29.*t62)-nl2.*(lambda.*t20.*t32.*t62+t10.*t20.*t32.*t62);zero;nl3.*(mu+lambda.*t20.*t70+t10.*t20.*t70)-nl2.*t66+nl1.*t69;nl2.*t73-nl1.*t76-nl3.*t79;zero;-nl3.*t79+nl2.*t82-nl1.*t85;nl3.*(mu+lambda.*t20.*t92+t10.*t20.*t92)-nl2.*t88+nl1.*t91;zero;nl1.*(q22.*t10.*t11+lambda.*t17.*t20.*t94+t10.*t17.*t20.*t94)-nl2.*(q21.*t10.*t11+lambda.*t20.*t22.*t94+t10.*t20.*t22.*t94)+nl3.*(lambda.*t15.*t20.*t94+t10.*t15.*t20.*t94);-nl1.*(q12.*t10.*t11+lambda.*t20.*t26.*t94+t10.*t20.*t26.*t94)+nl2.*(q11.*t10.*t11+lambda.*t20.*t32.*t94+t10.*t20.*t32.*t94)-nl3.*(lambda.*t20.*t29.*t94+t10.*t20.*t29.*t94);zero];
        end
        if nargout > 2
            fb_uh = [-t13;zero;zero;zero;-t13;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 6
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
        u3 = udg(:,3);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x1 = pg(:,1);
        x2 = pg(:,2);
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
        t13 = q12.*q23;
        t15 = q13.*q22;
        t14 = t13-t15;
        t16 = q22.*q33;
        t21 = q23.*q32;
        t17 = t16-t21;
        t18 = 1.0./(t2+t3+t4-t6-t7-t8).^2;
        t19 = q11.*q22;
        t22 = q12.*q21;
        t20 = t19-t22;
        t23 = q11.*q23;
        t25 = q13.*q21;
        t24 = t23-t25;
        fb = [one.*(-uh1+uinf1+x1);-uh2+uinf2+x2;-uinf3+nl1.*(mu.*q31+t10.*t11.*t14)+nl3.*(mu.*q33+t10.*t11.*t20)+nl2.*(mu.*q32-t10.*t11.*t24)+tau.*(u3-uh3)];
        if nargout > 1
            t26 = q12.*q33;
            t28 = q13.*q32;
            t27 = t26-t28;
            t29 = t14.^2;
            t30 = q21.*q33;
            t32 = q23.*q31;
            t31 = t30-t32;
            t33 = q23.*t10.*t11;
            t34 = q11.*q33;
            t36 = q13.*q31;
            t35 = t34-t36;
            t37 = q13.*t10.*t11;
            t38 = lambda.*t14.*t18.*t24;
            t39 = t10.*t14.*t18.*t24;
            t40 = t38+t39;
            t41 = t24.^2;
            t42 = q21.*q32;
            t45 = q22.*q31;
            t43 = t42-t45;
            t44 = q21.*t10.*t11;
            t46 = q22.*t10.*t11;
            t47 = q11.*q32;
            t50 = q12.*q31;
            t48 = t47-t50;
            t49 = q11.*t10.*t11;
            t51 = q12.*t10.*t11;
            t52 = lambda.*t18.*t20.*t24;
            t53 = t10.*t18.*t20.*t24;
            t54 = t52+t53;
            t55 = lambda.*t14.*t18.*t20;
            t56 = t10.*t14.*t18.*t20;
            t57 = t55+t56;
            t58 = t20.^2;
            t59 = one.*tau;
            fb_udg = [zero;zero;zero;zero;zero;zero;zero;zero;t59;zero;zero;nl1.*(lambda.*t14.*t17.*t18+t10.*t14.*t17.*t18)-nl2.*(t33+lambda.*t17.*t18.*t24+t10.*t17.*t18.*t24)+nl3.*(t46+lambda.*t17.*t18.*t20+t10.*t17.*t18.*t20);zero;zero;-nl1.*(lambda.*t14.*t18.*t27+t10.*t14.*t18.*t27)+nl2.*(t37+lambda.*t18.*t24.*t27+t10.*t18.*t24.*t27)-nl3.*(t51+lambda.*t18.*t20.*t27+t10.*t18.*t20.*t27);zero;zero;nl1.*(mu+lambda.*t18.*t29+t10.*t18.*t29)-nl2.*t40+nl3.*t57;zero;zero;-nl1.*(-t33+lambda.*t14.*t18.*t31+t10.*t14.*t18.*t31)+nl2.*(lambda.*t18.*t24.*t31+t10.*t18.*t24.*t31)-nl3.*(t44+lambda.*t18.*t20.*t31+t10.*t18.*t20.*t31);zero;zero;nl1.*(-t37+lambda.*t14.*t18.*t35+t10.*t14.*t18.*t35)-nl2.*(lambda.*t18.*t24.*t35+t10.*t18.*t24.*t35)+nl3.*(t49+lambda.*t18.*t20.*t35+t10.*t18.*t20.*t35);zero;zero;nl2.*(mu+lambda.*t18.*t41+t10.*t18.*t41)-nl1.*t40-nl3.*t54;zero;zero;nl1.*(-t46+lambda.*t14.*t18.*t43+t10.*t14.*t18.*t43)-nl2.*(-t44+lambda.*t18.*t24.*t43+t10.*t18.*t24.*t43)+nl3.*(lambda.*t18.*t20.*t43+t10.*t18.*t20.*t43);zero;zero;-nl1.*(-t51+lambda.*t14.*t18.*t48+t10.*t14.*t18.*t48)+nl2.*(-t49+lambda.*t18.*t24.*t48+t10.*t18.*t24.*t48)-nl3.*(lambda.*t18.*t20.*t48+t10.*t18.*t20.*t48);zero;zero;nl3.*(mu+lambda.*t18.*t58+t10.*t18.*t58)-nl2.*t54+nl1.*t57];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-one;zero;zero;zero;-t59];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 7
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
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x2 = pg(:,2);
        x3 = pg(:,3);
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
        t13 = q21.*q32;
        t17 = q22.*q31;
        t14 = t13-t17;
        t15 = q22.*q33;
        t18 = q23.*q32;
        t16 = t15-t18;
        t19 = 1.0./(t2+t3+t4-t6-t7-t8).^2;
        t20 = q21.*q33;
        t22 = q23.*q31;
        t21 = t20-t22;
        fb = [-uinf1+nl1.*(mu.*q11+t10.*t11.*t16)+nl3.*(mu.*q13+t10.*t11.*t14)+nl2.*(mu.*q12-t10.*t11.*t21)+tau.*(u1-uh1);-uh2+uinf2+x2;-uh3+uinf3+x3];
        if nargout > 1
            t23 = t16.^2;
            t24 = q12.*q33;
            t26 = q13.*q32;
            t25 = t24-t26;
            t27 = q12.*q23;
            t29 = q13.*q22;
            t28 = t27-t29;
            t30 = lambda.*t16.*t19.*t21;
            t31 = t10.*t16.*t19.*t21;
            t32 = t30+t31;
            t33 = t21.^2;
            t34 = q11.*q33;
            t36 = q13.*q31;
            t35 = t34-t36;
            t37 = q11.*q23;
            t39 = q13.*q21;
            t38 = t37-t39;
            t40 = lambda.*t14.*t19.*t21;
            t41 = t10.*t14.*t19.*t21;
            t42 = t40+t41;
            t43 = lambda.*t14.*t16.*t19;
            t44 = t10.*t14.*t16.*t19;
            t45 = t43+t44;
            t46 = t14.^2;
            t47 = q11.*q32;
            t49 = q12.*q31;
            t48 = t47-t49;
            t50 = q11.*q22;
            t52 = q12.*q21;
            t51 = t50-t52;
            t53 = one.*tau;
            fb_udg = [t53;zero;zero;zero;zero;zero;zero;zero;zero;nl1.*(mu+lambda.*t19.*t23+t10.*t19.*t23)-nl2.*t32+nl3.*t45;zero;zero;-nl3.*(-q32.*t10.*t11+lambda.*t14.*t19.*t25+t10.*t14.*t19.*t25)+nl2.*(-q33.*t10.*t11+lambda.*t19.*t21.*t25+t10.*t19.*t21.*t25)-nl1.*(lambda.*t16.*t19.*t25+t10.*t16.*t19.*t25);zero;zero;nl3.*(-q22.*t10.*t11+lambda.*t14.*t19.*t28+t10.*t14.*t19.*t28)-nl2.*(-q23.*t10.*t11+lambda.*t19.*t21.*t28+t10.*t19.*t21.*t28)+nl1.*(lambda.*t16.*t19.*t28+t10.*t16.*t19.*t28);zero;zero;nl2.*(mu+lambda.*t19.*t33+t10.*t19.*t33)-nl1.*t32-nl3.*t42;zero;zero;nl3.*(-q31.*t10.*t11+lambda.*t14.*t19.*t35+t10.*t14.*t19.*t35)+nl1.*(q33.*t10.*t11+lambda.*t16.*t19.*t35+t10.*t16.*t19.*t35)-nl2.*(lambda.*t19.*t21.*t35+t10.*t19.*t21.*t35);zero;zero;-nl3.*(-q21.*t10.*t11+lambda.*t14.*t19.*t38+t10.*t14.*t19.*t38)-nl1.*(q23.*t10.*t11+lambda.*t16.*t19.*t38+t10.*t16.*t19.*t38)+nl2.*(lambda.*t19.*t21.*t38+t10.*t19.*t21.*t38);zero;zero;nl3.*(mu+lambda.*t19.*t46+t10.*t19.*t46)-nl2.*t42+nl1.*t45;zero;zero;-nl1.*(q32.*t10.*t11+lambda.*t16.*t19.*t48+t10.*t16.*t19.*t48)+nl2.*(q31.*t10.*t11+lambda.*t19.*t21.*t48+t10.*t19.*t21.*t48)-nl3.*(lambda.*t14.*t19.*t48+t10.*t14.*t19.*t48);zero;zero;nl1.*(q22.*t10.*t11+lambda.*t16.*t19.*t51+t10.*t16.*t19.*t51)-nl2.*(q21.*t10.*t11+lambda.*t19.*t21.*t51+t10.*t19.*t21.*t51)+nl3.*(lambda.*t14.*t19.*t51+t10.*t14.*t19.*t51);zero;zero];
        end
        if nargout > 2
            fb_uh = [-t53;zero;zero;zero;-one;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 8
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
        u2 = udg(:,2);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uh3 = uh(:,3);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        uinf3 = uinf(:,3);
        x1 = pg(:,1);
        x3 = pg(:,3);
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
        t13 = q12.*q33;
        t15 = q13.*q32;
        t14 = t13-t15;
        t16 = q22.*q33;
        t21 = q23.*q32;
        t17 = t16-t21;
        t18 = 1.0./(t2+t3+t4-t6-t7-t8).^2;
        t19 = q11.*q32;
        t22 = q12.*q31;
        t20 = t19-t22;
        t23 = q11.*q33;
        t25 = q13.*q31;
        t24 = t23-t25;
        fb = [one.*(-uh1+uinf1+x1);-uinf2+nl1.*(mu.*q21-t10.*t11.*t14)+nl3.*(mu.*q23-t10.*t11.*t20)+nl2.*(mu.*q22+t10.*t11.*t24)+tau.*(u2-uh2);-uh3+uinf3+x3];
        if nargout > 1
            t26 = t14.^2;
            t27 = q12.*q23;
            t29 = q13.*q22;
            t28 = t27-t29;
            t30 = q21.*q33;
            t32 = q23.*q31;
            t31 = t30-t32;
            t33 = q33.*t10.*t11;
            t34 = lambda.*t14.*t18.*t24;
            t35 = t10.*t14.*t18.*t24;
            t36 = t34+t35;
            t37 = t24.^2;
            t38 = q11.*q23;
            t40 = q13.*q21;
            t39 = t38-t40;
            t41 = q21.*q32;
            t44 = q22.*q31;
            t42 = t41-t44;
            t43 = q31.*t10.*t11;
            t45 = q32.*t10.*t11;
            t46 = lambda.*t18.*t20.*t24;
            t47 = t10.*t18.*t20.*t24;
            t48 = t46+t47;
            t49 = lambda.*t14.*t18.*t20;
            t50 = t10.*t14.*t18.*t20;
            t51 = t49+t50;
            t52 = t20.^2;
            t53 = q11.*q22;
            t55 = q12.*q21;
            t54 = t53-t55;
            t56 = one.*tau;
            fb_udg = [zero;zero;zero;zero;t56;zero;zero;zero;zero;zero;-nl1.*(lambda.*t14.*t17.*t18+t10.*t14.*t17.*t18)+nl2.*(t33+lambda.*t17.*t18.*t24+t10.*t17.*t18.*t24)-nl3.*(t45+lambda.*t17.*t18.*t20+t10.*t17.*t18.*t20);zero;zero;nl1.*(mu+lambda.*t18.*t26+t10.*t18.*t26)-nl2.*t36+nl3.*t51;zero;zero;-nl3.*(-q12.*t10.*t11+lambda.*t18.*t20.*t28+t10.*t18.*t20.*t28)+nl2.*(-q13.*t10.*t11+lambda.*t18.*t24.*t28+t10.*t18.*t24.*t28)-nl1.*(lambda.*t14.*t18.*t28+t10.*t14.*t18.*t28);zero;zero;nl1.*(-t33+lambda.*t14.*t18.*t31+t10.*t14.*t18.*t31)-nl2.*(lambda.*t18.*t24.*t31+t10.*t18.*t24.*t31)+nl3.*(t43+lambda.*t18.*t20.*t31+t10.*t18.*t20.*t31);zero;zero;nl2.*(mu+lambda.*t18.*t37+t10.*t18.*t37)-nl1.*t36-nl3.*t48;zero;zero;nl1.*(q13.*t10.*t11+lambda.*t14.*t18.*t39+t10.*t14.*t18.*t39)+nl3.*(-q11.*t10.*t11+lambda.*t18.*t20.*t39+t10.*t18.*t20.*t39)-nl2.*(lambda.*t18.*t24.*t39+t10.*t18.*t24.*t39);zero;zero;-nl1.*(-t45+lambda.*t14.*t18.*t42+t10.*t14.*t18.*t42)+nl2.*(-t43+lambda.*t18.*t24.*t42+t10.*t18.*t24.*t42)-nl3.*(lambda.*t18.*t20.*t42+t10.*t18.*t20.*t42);zero;zero;nl3.*(mu+lambda.*t18.*t52+t10.*t18.*t52)-nl2.*t48+nl1.*t51;zero;zero;-nl1.*(q12.*t10.*t11+lambda.*t14.*t18.*t54+t10.*t14.*t18.*t54)+nl2.*(q11.*t10.*t11+lambda.*t18.*t24.*t54+t10.*t18.*t24.*t54)-nl3.*(lambda.*t18.*t20.*t54+t10.*t18.*t20.*t54);zero];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-t56;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    otherwise
         error('unknown boundary type');
end
