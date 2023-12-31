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
        [fb,fb_udg,fb_uh] = fhat(nl,pg,udg,uh,param,time);
        fb = fb - uinf;
    case 3
        kappa = param{2};
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
        t2 = q11.*q32;
        t24 = q12.*q31;
        t3 = t2-t24;
        t4 = q11.*q23.*q32;
        t5 = q12.*q21.*q33;
        t6 = q13.*q22.*q31;
        t10 = q11.*q22.*q33;
        t11 = q12.*q23.*q31;
        t12 = q13.*q21.*q32;
        t7 = t4+t5+t6-t10-t11-t12+1.0;
        t8 = q11.*q33;
        t27 = q13.*q31;
        t9 = t8-t27;
        t13 = q12.*q33;
        t30 = q13.*q32;
        t14 = t13-t30;
        t15 = t4+t5+t6-t10-t11-t12;
        t16 = 1.0./t15;
        t17 = q11.*q22;
        t31 = q12.*q21;
        t18 = t17-t31;
        t19 = q11.*q23;
        t32 = q13.*q21;
        t20 = t19-t32;
        t21 = q12.*q23;
        t33 = q13.*q22;
        t22 = t21-t33;
        fb = [one.*(-uh1+uinf1+x1);-uinf2+tau.*(u2-uh2)+nl3.*(mu.*q23-(mu.*t3)./(t4+t5+t6-q11.*q22.*q33-q12.*q23.*q31-q13.*q21.*q32)+kappa.*t3.*t7)+nl2.*(mu.*q22-kappa.*t7.*t9+mu.*t9.*t16)+nl1.*(mu.*q21+kappa.*t7.*t14-mu.*t14.*t16);-uinf3+tau.*(u3-uh3)+nl3.*(mu.*q33-kappa.*t7.*t18+mu.*t16.*t18)+nl2.*(mu.*q32+kappa.*t7.*t20-mu.*t16.*t20)+nl1.*(mu.*q31-kappa.*t7.*t22+mu.*t16.*t22)];
        if nargout > 1
            t23 = one.*tau;
            t25 = q22.*q33;
            t28 = q23.*q32;
            t26 = t25-t28;
            t29 = 1.0./t15.^2;
            t34 = 1.0./(t4+t5+t6-t10-t11-t12).^2;
            t35 = t14.^2;
            t36 = mu.*q12.*t16;
            t37 = mu.*q13.*t16;
            t38 = kappa.*t14.*t22;
            t39 = mu.*t14.*t22.*t34;
            t40 = t38+t39;
            t41 = t22.^2;
            t42 = q21.*q33;
            t45 = q23.*q31;
            t43 = t42-t45;
            t44 = mu.*q33.*t16;
            t46 = mu.*q23.*t16;
            t47 = kappa.*t9.*t14;
            t48 = mu.*t9.*t14.*t34;
            t49 = t47+t48;
            t50 = t9.^2;
            t51 = kappa.*t9.*t22;
            t52 = kappa.*q13.*t7;
            t53 = mu.*t9.*t22.*t34;
            t54 = -t37+t51+t52+t53;
            t55 = mu.*q11.*t16;
            t56 = kappa.*t14.*t20;
            t57 = mu.*t14.*t20.*t34;
            t58 = kappa.*t9.*t20;
            t59 = mu.*t9.*t20.*t34;
            t60 = t58+t59;
            t61 = kappa.*t20.*t22;
            t62 = mu.*t20.*t22.*t34;
            t63 = t61+t62;
            t64 = t20.^2;
            t65 = mu.*q31.*t16;
            t66 = q21.*q32;
            t69 = q22.*q31;
            t67 = t66-t69;
            t68 = mu.*q32.*t16;
            t70 = mu.*q21.*t16;
            t71 = mu.*q22.*t16;
            t72 = kappa.*t3.*t9;
            t73 = mu.*t3.*t9.*t34;
            t74 = t72+t73;
            t75 = kappa.*t3.*t14;
            t76 = mu.*t3.*t14.*t34;
            t77 = t75+t76;
            t78 = t3.^2;
            t79 = kappa.*t3.*t20;
            t80 = kappa.*q11.*t7;
            t81 = mu.*t3.*t20.*t34;
            t82 = -t55+t79+t80+t81;
            t83 = kappa.*t3.*t22;
            t84 = kappa.*q12.*t7;
            t85 = mu.*t3.*t22.*t34;
            t86 = -t36+t83+t84+t85;
            t87 = kappa.*t9.*t18;
            t88 = mu.*t9.*t18.*t34;
            t89 = kappa.*t14.*t18;
            t90 = mu.*t14.*t18.*t34;
            t91 = kappa.*t3.*t18;
            t92 = mu.*t3.*t18.*t34;
            t93 = t91+t92;
            t94 = kappa.*t18.*t20;
            t95 = mu.*t18.*t20.*t34;
            t96 = t94+t95;
            t97 = kappa.*t18.*t22;
            t98 = mu.*t18.*t22.*t34;
            t99 = t97+t98;
            t100 = t18.^2;
            fb_udg = [zero;zero;zero;zero;t23;zero;zero;zero;t23;zero;-nl1.*(kappa.*t14.*t26+mu.*t14.*t26.*t29)+nl2.*(t44-kappa.*q33.*t7+kappa.*t9.*t26+mu.*t9.*t26.*t29)-nl3.*(t68-kappa.*q32.*t7+kappa.*t3.*t26+mu.*t3.*t26.*t29);nl1.*(kappa.*t22.*t26+mu.*t22.*t26.*t34)-nl2.*(t46-kappa.*q23.*t7+kappa.*t20.*t26+mu.*t20.*t26.*t29)+nl3.*(t71-kappa.*q22.*t7+kappa.*t18.*t26+mu.*t18.*t26.*t29);zero;-nl2.*t49+nl3.*t77+nl1.*(mu+kappa.*t35+mu.*t34.*t35);nl2.*(t37+t56+t57-kappa.*q13.*t7)-nl3.*(t36+t89+t90-kappa.*q12.*t7)-nl1.*t40;zero;-nl1.*t40+nl2.*t54-nl3.*t86;-nl2.*t63+nl3.*t99+nl1.*(mu+kappa.*t41+mu.*t34.*t41);zero;-nl2.*(kappa.*t9.*t43+mu.*t9.*t34.*t43)+nl3.*(t65-kappa.*q31.*t7+kappa.*t3.*t43+mu.*t3.*t34.*t43)+nl1.*(-t44+kappa.*q33.*t7+kappa.*t14.*t43+mu.*t14.*t34.*t43);nl2.*(kappa.*t20.*t43+mu.*t20.*t34.*t43)-nl3.*(t70-kappa.*q21.*t7+kappa.*t18.*t43+mu.*t18.*t34.*t43)-nl1.*(-t46+kappa.*q23.*t7+kappa.*t22.*t43+mu.*t22.*t34.*t43);zero;-nl1.*t49-nl3.*t74+nl2.*(mu+kappa.*t50+mu.*t34.*t50);nl3.*(t55+t87+t88-kappa.*q11.*t7)+nl1.*t54-nl2.*t60;zero;-nl2.*t60+nl3.*t82+nl1.*(t37-t52+t56+t57);-nl1.*t63-nl3.*t96+nl2.*(mu+kappa.*t64+mu.*t34.*t64);zero;-nl3.*(kappa.*t3.*t67+mu.*t3.*t34.*t67)+nl2.*(-t65+kappa.*q31.*t7+kappa.*t9.*t67+mu.*t9.*t34.*t67)-nl1.*(-t68+kappa.*q32.*t7+kappa.*t14.*t67+mu.*t14.*t34.*t67);nl3.*(kappa.*t18.*t67+mu.*t18.*t34.*t67)-nl2.*(-t70+kappa.*q21.*t7+kappa.*t20.*t67+mu.*t20.*t34.*t67)+nl1.*(-t71+kappa.*q22.*t7+kappa.*t22.*t67+mu.*t22.*t34.*t67);zero;-nl2.*t74+nl1.*t77+nl3.*(mu+kappa.*t78+mu.*t34.*t78);nl2.*t82-nl1.*t86-nl3.*t93;zero;-nl3.*t93-nl1.*(t36-t84+t89+t90)+nl2.*(t55-t80+t87+t88);-nl2.*t96+nl1.*t99+nl3.*(mu+kappa.*t100+mu.*t34.*t100)];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-t23;zero;zero;zero;-t23];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 4
        kappa = param{2};
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
        t2 = q21.*q32;
        t24 = q22.*q31;
        t3 = t2-t24;
        t4 = q11.*q23.*q32;
        t5 = q12.*q21.*q33;
        t6 = q13.*q22.*q31;
        t10 = q11.*q22.*q33;
        t11 = q12.*q23.*q31;
        t12 = q13.*q21.*q32;
        t7 = t4+t5+t6-t10-t11-t12+1.0;
        t8 = q21.*q33;
        t26 = q23.*q31;
        t9 = t8-t26;
        t13 = q22.*q33;
        t25 = q23.*q32;
        t14 = t13-t25;
        t15 = t4+t5+t6-t10-t11-t12;
        t16 = 1.0./t15;
        t17 = q11.*q22;
        t29 = q12.*q21;
        t18 = t17-t29;
        t19 = q11.*q23;
        t30 = q13.*q21;
        t20 = t19-t30;
        t21 = q12.*q23;
        t31 = q13.*q22;
        t22 = t21-t31;
        fb = [-uinf1+tau.*(u1-uh1)+nl3.*(mu.*q13+(mu.*t3)./(t4+t5+t6-q11.*q22.*q33-q12.*q23.*q31-q13.*q21.*q32)-kappa.*t3.*t7)+nl2.*(mu.*q12+kappa.*t7.*t9-mu.*t9.*t16)+nl1.*(mu.*q11-kappa.*t7.*t14+mu.*t14.*t16);-uh2+uinf2+x2;-uinf3+tau.*(u3-uh3)+nl3.*(mu.*q33-kappa.*t7.*t18+mu.*t16.*t18)+nl2.*(mu.*q32+kappa.*t7.*t20-mu.*t16.*t20)+nl1.*(mu.*q31-kappa.*t7.*t22+mu.*t16.*t22)];
        if nargout > 1
            t23 = one.*tau;
            t27 = 1.0./(t4+t5+t6-t10-t11-t12).^2;
            t28 = t14.^2;
            t32 = q12.*q33;
            t34 = q13.*q32;
            t33 = t32-t34;
            t35 = mu.*q22.*t16;
            t36 = mu.*q23.*t16;
            t37 = kappa.*t14.*t22;
            t38 = mu.*t14.*t22.*t27;
            t39 = t37+t38;
            t40 = nl1.*t39;
            t41 = t22.^2;
            t42 = kappa.*t9.*t14;
            t43 = mu.*t9.*t14.*t27;
            t44 = t42+t43;
            t45 = t9.^2;
            t46 = kappa.*t9.*t22;
            t47 = kappa.*q23.*t7;
            t48 = mu.*t9.*t22.*t27;
            t49 = -t36+t46+t47+t48;
            t50 = q11.*q33;
            t52 = q13.*q31;
            t51 = t50-t52;
            t53 = kappa.*q33.*t7;
            t54 = mu.*q13.*t16;
            t55 = mu.*q21.*t16;
            t56 = kappa.*t14.*t20;
            t57 = mu.*t14.*t20.*t27;
            t58 = kappa.*t9.*t20;
            t59 = mu.*t9.*t20.*t27;
            t60 = t58+t59;
            t61 = nl2.*t60;
            t62 = kappa.*t20.*t22;
            t63 = mu.*t20.*t22.*t27;
            t64 = t62+t63;
            t65 = t20.^2;
            t66 = kappa.*t3.*t9;
            t67 = mu.*t3.*t9.*t27;
            t68 = t66+t67;
            t69 = kappa.*t3.*t14;
            t70 = mu.*t3.*t14.*t27;
            t71 = t69+t70;
            t72 = t3.^2;
            t73 = kappa.*t3.*t20;
            t74 = kappa.*q21.*t7;
            t75 = mu.*t3.*t20.*t27;
            t76 = -t55+t73+t74+t75;
            t77 = kappa.*t3.*t22;
            t78 = kappa.*q22.*t7;
            t79 = mu.*t3.*t22.*t27;
            t80 = -t35+t77+t78+t79;
            t81 = kappa.*q31.*t7;
            t82 = q11.*q32;
            t84 = q12.*q31;
            t83 = t82-t84;
            t85 = kappa.*q32.*t7;
            t86 = mu.*q11.*t16;
            t87 = mu.*q12.*t16;
            t88 = kappa.*t9.*t18;
            t89 = mu.*t9.*t18.*t27;
            t90 = kappa.*t14.*t18;
            t91 = mu.*t14.*t18.*t27;
            t92 = kappa.*t3.*t18;
            t93 = mu.*t3.*t18.*t27;
            t94 = t92+t93;
            t95 = nl3.*t94;
            t96 = kappa.*t18.*t20;
            t97 = mu.*t18.*t20.*t27;
            t98 = t96+t97;
            t99 = kappa.*t18.*t22;
            t100 = mu.*t18.*t22.*t27;
            t101 = t99+t100;
            t102 = t18.^2;
            fb_udg = [t23;zero;zero;zero;zero;zero;zero;zero;t23;-nl2.*t44+nl3.*t71+nl1.*(mu+kappa.*t28+mu.*t27.*t28);zero;t40-nl2.*(t36+t56+t57-kappa.*q23.*t7)+nl3.*(t35+t90+t91-kappa.*q22.*t7);-nl1.*(kappa.*t14.*t33+mu.*t14.*t27.*t33)+nl2.*(t53-mu.*q33.*t16+kappa.*t9.*t33+mu.*t9.*t27.*t33)-nl3.*(t85-mu.*q32.*t16+kappa.*t3.*t33+mu.*t3.*t27.*t33);zero;-nl1.*(kappa.*t22.*t33+mu.*t22.*t27.*t33)+nl2.*(t54-kappa.*q13.*t7+kappa.*t20.*t33+mu.*t20.*t27.*t33)-nl3.*(t87-kappa.*q12.*t7+kappa.*t18.*t33+mu.*t18.*t27.*t33);t40-nl2.*t49+nl3.*t80;zero;-nl2.*t64+nl3.*t101+nl1.*(mu+kappa.*t41+mu.*t27.*t41);-nl1.*t44-nl3.*t68+nl2.*(mu+kappa.*t45+mu.*t27.*t45);zero;t61-nl3.*(t55+t88+t89-kappa.*q21.*t7)-nl1.*t49;-nl2.*(kappa.*t9.*t51+mu.*t9.*t27.*t51)+nl3.*(t81-mu.*q31.*t16+kappa.*t3.*t51+mu.*t3.*t27.*t51)+nl1.*(-t53+mu.*q33.*t16+kappa.*t14.*t51+mu.*t14.*t27.*t51);zero;-nl2.*(kappa.*t20.*t51+mu.*t20.*t27.*t51)+nl3.*(t86-kappa.*q11.*t7+kappa.*t18.*t51+mu.*t18.*t27.*t51)+nl1.*(-t54+kappa.*q13.*t7+kappa.*t22.*t51+mu.*t22.*t27.*t51);t61-nl3.*t76-nl1.*(t36-t47+t56+t57);zero;-nl1.*t64-nl3.*t98+nl2.*(mu+kappa.*t65+mu.*t27.*t65);-nl2.*t68+nl1.*t71+nl3.*(mu+kappa.*t72+mu.*t27.*t72);zero;t95-nl2.*t76+nl1.*t80;-nl3.*(kappa.*t3.*t83+mu.*t3.*t27.*t83)+nl2.*(-t81+mu.*q31.*t16+kappa.*t9.*t83+mu.*t9.*t27.*t83)-nl1.*(-t85+mu.*q32.*t16+kappa.*t14.*t83+mu.*t14.*t27.*t83);zero;-nl3.*(kappa.*t18.*t83+mu.*t18.*t27.*t83)+nl2.*(-t86+kappa.*q11.*t7+kappa.*t20.*t83+mu.*t20.*t27.*t83)-nl1.*(-t87+kappa.*q12.*t7+kappa.*t22.*t83+mu.*t22.*t27.*t83);t95+nl1.*(t35-t78+t90+t91)-nl2.*(t55-t74+t88+t89);zero;-nl2.*t98+nl1.*t101+nl3.*(mu+kappa.*t102+mu.*t27.*t102)];
        end
        if nargout > 2
            fb_uh = [-t23;zero;zero;zero;-one;zero;zero;zero;-t23];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 5
        kappa = param{2};
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
        t2 = q21.*q32;
        t24 = q22.*q31;
        t3 = t2-t24;
        t4 = q11.*q23.*q32;
        t5 = q12.*q21.*q33;
        t6 = q13.*q22.*q31;
        t10 = q11.*q22.*q33;
        t11 = q12.*q23.*q31;
        t12 = q13.*q21.*q32;
        t7 = t4+t5+t6-t10-t11-t12+1.0;
        t8 = q21.*q33;
        t26 = q23.*q31;
        t9 = t8-t26;
        t13 = q22.*q33;
        t25 = q23.*q32;
        t14 = t13-t25;
        t15 = t4+t5+t6-t10-t11-t12;
        t16 = 1.0./t15;
        t17 = q11.*q32;
        t29 = q12.*q31;
        t18 = t17-t29;
        t19 = q11.*q33;
        t30 = q13.*q31;
        t20 = t19-t30;
        t21 = q12.*q33;
        t31 = q13.*q32;
        t22 = t21-t31;
        fb = [-uinf1+tau.*(u1-uh1)+nl3.*(mu.*q13+(mu.*t3)./(t4+t5+t6-q11.*q22.*q33-q12.*q23.*q31-q13.*q21.*q32)-kappa.*t3.*t7)+nl2.*(mu.*q12+kappa.*t7.*t9-mu.*t9.*t16)+nl1.*(mu.*q11-kappa.*t7.*t14+mu.*t14.*t16);-uinf2+tau.*(u2-uh2)+nl3.*(mu.*q23+kappa.*t7.*t18-mu.*t16.*t18)+nl2.*(mu.*q22-kappa.*t7.*t20+mu.*t16.*t20)+nl1.*(mu.*q21+kappa.*t7.*t22-mu.*t16.*t22);-uh3+uinf3+x3];
        if nargout > 1
            t23 = one.*tau;
            t27 = 1.0./(t4+t5+t6-t10-t11-t12).^2;
            t28 = t14.^2;
            t32 = mu.*q32.*t16;
            t33 = mu.*q33.*t16;
            t34 = kappa.*t14.*t22;
            t35 = mu.*t14.*t22.*t27;
            t36 = t34+t35;
            t37 = t22.^2;
            t38 = q12.*q23;
            t40 = q13.*q22;
            t39 = t38-t40;
            t41 = kappa.*t9.*t14;
            t42 = mu.*t9.*t14.*t27;
            t43 = t41+t42;
            t44 = t9.^2;
            t45 = kappa.*t9.*t22;
            t46 = kappa.*q33.*t7;
            t47 = mu.*t9.*t22.*t27;
            t48 = -t33+t45+t46+t47;
            t49 = mu.*q31.*t16;
            t50 = kappa.*t14.*t20;
            t51 = mu.*t14.*t20.*t27;
            t52 = kappa.*t9.*t20;
            t53 = mu.*t9.*t20.*t27;
            t54 = t52+t53;
            t55 = kappa.*t20.*t22;
            t56 = mu.*t20.*t22.*t27;
            t57 = t55+t56;
            t58 = t20.^2;
            t59 = q11.*q23;
            t61 = q13.*q21;
            t60 = t59-t61;
            t62 = kappa.*q23.*t7;
            t63 = kappa.*q13.*t7;
            t64 = kappa.*t3.*t9;
            t65 = mu.*t3.*t9.*t27;
            t66 = t64+t65;
            t67 = kappa.*t3.*t14;
            t68 = mu.*t3.*t14.*t27;
            t69 = t67+t68;
            t70 = t3.^2;
            t71 = kappa.*t3.*t20;
            t72 = kappa.*q31.*t7;
            t73 = mu.*t3.*t20.*t27;
            t74 = -t49+t71+t72+t73;
            t75 = kappa.*t3.*t22;
            t76 = kappa.*q32.*t7;
            t77 = mu.*t3.*t22.*t27;
            t78 = -t32+t75+t76+t77;
            t79 = kappa.*t9.*t18;
            t80 = mu.*t9.*t18.*t27;
            t81 = kappa.*t14.*t18;
            t82 = mu.*t14.*t18.*t27;
            t83 = kappa.*t3.*t18;
            t84 = mu.*t3.*t18.*t27;
            t85 = t83+t84;
            t86 = kappa.*t18.*t20;
            t87 = mu.*t18.*t20.*t27;
            t88 = t86+t87;
            t89 = kappa.*t18.*t22;
            t90 = mu.*t18.*t22.*t27;
            t91 = t89+t90;
            t92 = t18.^2;
            t93 = kappa.*q21.*t7;
            t94 = q11.*q22;
            t96 = q12.*q21;
            t95 = t94-t96;
            t97 = kappa.*q22.*t7;
            t98 = kappa.*q11.*t7;
            t99 = kappa.*q12.*t7;
            fb_udg = [t23;zero;zero;zero;t23;zero;zero;zero;zero;-nl2.*t43+nl3.*t69+nl1.*(mu+kappa.*t28+mu.*t27.*t28);nl2.*(t33+t50+t51-kappa.*q33.*t7)-nl3.*(t32+t81+t82-kappa.*q32.*t7)-nl1.*t36;zero;-nl1.*t36+nl2.*t48-nl3.*t78;-nl2.*t57+nl3.*t91+nl1.*(mu+kappa.*t37+mu.*t27.*t37);zero;nl1.*(kappa.*t14.*t39+mu.*t14.*t27.*t39)-nl2.*(t62-mu.*q23.*t16+kappa.*t9.*t39+mu.*t9.*t27.*t39)+nl3.*(t97-mu.*q22.*t16+kappa.*t3.*t39+mu.*t3.*t27.*t39);-nl1.*(kappa.*t22.*t39+mu.*t22.*t27.*t39)+nl2.*(t63-mu.*q13.*t16+kappa.*t20.*t39+mu.*t20.*t27.*t39)-nl3.*(t99-mu.*q12.*t16+kappa.*t18.*t39+mu.*t18.*t27.*t39);zero;-nl1.*t43-nl3.*t66+nl2.*(mu+kappa.*t44+mu.*t27.*t44);nl3.*(t49+t79+t80-kappa.*q31.*t7)+nl1.*t48-nl2.*t54;zero;-nl2.*t54+nl3.*t74+nl1.*(t33-t46+t50+t51);-nl1.*t57-nl3.*t88+nl2.*(mu+kappa.*t58+mu.*t27.*t58);zero;nl2.*(kappa.*t9.*t60+mu.*t9.*t27.*t60)-nl3.*(t93-mu.*q21.*t16+kappa.*t3.*t60+mu.*t3.*t27.*t60)-nl1.*(-t62+mu.*q23.*t16+kappa.*t14.*t60+mu.*t14.*t27.*t60);-nl2.*(kappa.*t20.*t60+mu.*t20.*t27.*t60)+nl3.*(t98-mu.*q11.*t16+kappa.*t18.*t60+mu.*t18.*t27.*t60)+nl1.*(-t63+mu.*q13.*t16+kappa.*t22.*t60+mu.*t22.*t27.*t60);zero;-nl2.*t66+nl1.*t69+nl3.*(mu+kappa.*t70+mu.*t27.*t70);nl2.*t74-nl1.*t78-nl3.*t85;zero;-nl3.*t85-nl1.*(t32-t76+t81+t82)+nl2.*(t49-t72+t79+t80);-nl2.*t88+nl1.*t91+nl3.*(mu+kappa.*t92+mu.*t27.*t92);zero;nl3.*(kappa.*t3.*t95+mu.*t3.*t27.*t95)-nl2.*(-t93+mu.*q21.*t16+kappa.*t9.*t95+mu.*t9.*t27.*t95)+nl1.*(-t97+mu.*q22.*t16+kappa.*t14.*t95+mu.*t14.*t27.*t95);-nl3.*(kappa.*t18.*t95+mu.*t18.*t27.*t95)+nl2.*(-t98+mu.*q11.*t16+kappa.*t20.*t95+mu.*t20.*t27.*t95)-nl1.*(-t99+mu.*q12.*t16+kappa.*t22.*t95+mu.*t22.*t27.*t95);zero];
        end
        if nargout > 2
            fb_uh = [-t23;zero;zero;zero;-t23;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 6
        kappa = param{2};
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
        t2 = q11.*q22;
        t17 = q12.*q21;
        t3 = t2-t17;
        t4 = q11.*q23.*q32;
        t5 = q12.*q21.*q33;
        t6 = q13.*q22.*q31;
        t10 = q11.*q22.*q33;
        t11 = q12.*q23.*q31;
        t12 = q13.*q21.*q32;
        t7 = t4+t5+t6-t10-t11-t12+1.0;
        t8 = q11.*q23;
        t20 = q13.*q21;
        t9 = t8-t20;
        t13 = q12.*q23;
        t23 = q13.*q22;
        t14 = t13-t23;
        t15 = t4+t5+t6-t10-t11-t12;
        t16 = 1.0./t15;
        fb = [one.*(-uh1+uinf1+x1);-uh2+uinf2+x2;-uinf3+tau.*(u3-uh3)+nl3.*(mu.*q33+(mu.*t3)./(t4+t5+t6-q11.*q22.*q33-q12.*q23.*q31-q13.*q21.*q32)-kappa.*t3.*t7)+nl2.*(mu.*q32+kappa.*t7.*t9-mu.*t9.*t16)+nl1.*(mu.*q31-kappa.*t7.*t14+mu.*t14.*t16)];
        if nargout > 1
            t18 = q22.*q33;
            t21 = q23.*q32;
            t19 = t18-t21;
            t22 = 1.0./t15.^2;
            t24 = q12.*q33;
            t27 = q13.*q32;
            t25 = t24-t27;
            t26 = 1.0./(t4+t5+t6-t10-t11-t12).^2;
            t28 = t14.^2;
            t29 = q21.*q33;
            t32 = q23.*q31;
            t30 = t29-t32;
            t31 = mu.*q23.*t16;
            t33 = q11.*q33;
            t36 = q13.*q31;
            t34 = t33-t36;
            t35 = mu.*q13.*t16;
            t37 = kappa.*t9.*t14;
            t38 = mu.*t9.*t14.*t26;
            t39 = t37+t38;
            t40 = t9.^2;
            t41 = mu.*q21.*t16;
            t42 = q21.*q32;
            t45 = q22.*q31;
            t43 = t42-t45;
            t44 = mu.*q22.*t16;
            t46 = mu.*q11.*t16;
            t47 = q11.*q32;
            t50 = q12.*q31;
            t48 = t47-t50;
            t49 = mu.*q12.*t16;
            t51 = kappa.*t3.*t9;
            t52 = mu.*t3.*t9.*t26;
            t53 = t51+t52;
            t54 = kappa.*t3.*t14;
            t55 = mu.*t3.*t14.*t26;
            t56 = t54+t55;
            t57 = t3.^2;
            t58 = one.*tau;
            fb_udg = [zero;zero;zero;zero;zero;zero;zero;zero;t58;zero;zero;nl1.*(kappa.*t14.*t19+mu.*t14.*t19.*t26)-nl2.*(t31-kappa.*q23.*t7+kappa.*t9.*t19+mu.*t9.*t19.*t22)+nl3.*(t44-kappa.*q22.*t7+kappa.*t3.*t19+mu.*t3.*t19.*t22);zero;zero;-nl1.*(kappa.*t14.*t25+mu.*t14.*t25.*t26)+nl2.*(t35-kappa.*q13.*t7+kappa.*t9.*t25+mu.*t9.*t25.*t26)-nl3.*(t49-kappa.*q12.*t7+kappa.*t3.*t25+mu.*t3.*t25.*t26);zero;zero;-nl2.*t39+nl3.*t56+nl1.*(mu+kappa.*t28+mu.*t26.*t28);zero;zero;nl2.*(kappa.*t9.*t30+mu.*t9.*t26.*t30)-nl3.*(t41-kappa.*q21.*t7+kappa.*t3.*t30+mu.*t3.*t26.*t30)-nl1.*(-t31+kappa.*q23.*t7+kappa.*t14.*t30+mu.*t14.*t26.*t30);zero;zero;-nl2.*(kappa.*t9.*t34+mu.*t9.*t26.*t34)+nl3.*(t46-kappa.*q11.*t7+kappa.*t3.*t34+mu.*t3.*t26.*t34)+nl1.*(-t35+kappa.*q13.*t7+kappa.*t14.*t34+mu.*t14.*t26.*t34);zero;zero;-nl1.*t39-nl3.*t53+nl2.*(mu+kappa.*t40+mu.*t26.*t40);zero;zero;nl3.*(kappa.*t3.*t43+mu.*t3.*t26.*t43)-nl2.*(-t41+kappa.*q21.*t7+kappa.*t9.*t43+mu.*t9.*t26.*t43)+nl1.*(-t44+kappa.*q22.*t7+kappa.*t14.*t43+mu.*t14.*t26.*t43);zero;zero;-nl3.*(kappa.*t3.*t48+mu.*t3.*t26.*t48)+nl2.*(-t46+kappa.*q11.*t7+kappa.*t9.*t48+mu.*t9.*t26.*t48)-nl1.*(-t49+kappa.*q12.*t7+kappa.*t14.*t48+mu.*t14.*t26.*t48);zero;zero;-nl2.*t53+nl1.*t56+nl3.*(mu+kappa.*t57+mu.*t26.*t57)];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-one;zero;zero;zero;-t58];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 7
        kappa = param{2};
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
        t2 = q21.*q32;
        t17 = q22.*q31;
        t3 = t2-t17;
        t4 = q11.*q23.*q32;
        t5 = q12.*q21.*q33;
        t6 = q13.*q22.*q31;
        t10 = q11.*q22.*q33;
        t11 = q12.*q23.*q31;
        t12 = q13.*q21.*q32;
        t7 = t4+t5+t6-t10-t11-t12+1.0;
        t8 = q21.*q33;
        t19 = q23.*q31;
        t9 = t8-t19;
        t13 = q22.*q33;
        t18 = q23.*q32;
        t14 = t13-t18;
        t15 = t4+t5+t6-t10-t11-t12;
        t16 = 1.0./t15;
        fb = [-uinf1+tau.*(u1-uh1)+nl3.*(mu.*q13+(mu.*t3)./(t4+t5+t6-q11.*q22.*q33-q12.*q23.*q31-q13.*q21.*q32)-kappa.*t3.*t7)+nl2.*(mu.*q12+kappa.*t7.*t9-mu.*t9.*t16)+nl1.*(mu.*q11-kappa.*t7.*t14+mu.*t14.*t16);-uh2+uinf2+x2;-uh3+uinf3+x3];
        if nargout > 1
            t20 = 1.0./(t4+t5+t6-t10-t11-t12).^2;
            t21 = t14.^2;
            t22 = q12.*q33;
            t24 = q13.*q32;
            t23 = t22-t24;
            t25 = q12.*q23;
            t27 = q13.*q22;
            t26 = t25-t27;
            t28 = kappa.*t9.*t14;
            t29 = mu.*t9.*t14.*t20;
            t30 = t28+t29;
            t31 = t9.^2;
            t32 = q11.*q33;
            t34 = q13.*q31;
            t33 = t32-t34;
            t35 = kappa.*q33.*t7;
            t36 = q11.*q23;
            t38 = q13.*q21;
            t37 = t36-t38;
            t39 = kappa.*q23.*t7;
            t40 = kappa.*t3.*t9;
            t41 = mu.*t3.*t9.*t20;
            t42 = t40+t41;
            t43 = kappa.*t3.*t14;
            t44 = mu.*t3.*t14.*t20;
            t45 = t43+t44;
            t46 = t3.^2;
            t47 = kappa.*q31.*t7;
            t48 = q11.*q32;
            t50 = q12.*q31;
            t49 = t48-t50;
            t51 = kappa.*q32.*t7;
            t52 = kappa.*q21.*t7;
            t53 = q11.*q22;
            t55 = q12.*q21;
            t54 = t53-t55;
            t56 = kappa.*q22.*t7;
            t57 = one.*tau;
            fb_udg = [t57;zero;zero;zero;zero;zero;zero;zero;zero;-nl2.*t30+nl3.*t45+nl1.*(mu+kappa.*t21+mu.*t20.*t21);zero;zero;-nl1.*(kappa.*t14.*t23+mu.*t14.*t20.*t23)+nl2.*(t35-mu.*q33.*t16+kappa.*t9.*t23+mu.*t9.*t20.*t23)-nl3.*(t51-mu.*q32.*t16+kappa.*t3.*t23+mu.*t3.*t20.*t23);zero;zero;nl1.*(kappa.*t14.*t26+mu.*t14.*t20.*t26)-nl2.*(t39-mu.*q23.*t16+kappa.*t9.*t26+mu.*t9.*t20.*t26)+nl3.*(t56-mu.*q22.*t16+kappa.*t3.*t26+mu.*t3.*t20.*t26);zero;zero;-nl1.*t30-nl3.*t42+nl2.*(mu+kappa.*t31+mu.*t20.*t31);zero;zero;-nl2.*(kappa.*t9.*t33+mu.*t9.*t20.*t33)+nl3.*(t47-mu.*q31.*t16+kappa.*t3.*t33+mu.*t3.*t20.*t33)+nl1.*(-t35+mu.*q33.*t16+kappa.*t14.*t33+mu.*t14.*t20.*t33);zero;zero;nl2.*(kappa.*t9.*t37+mu.*t9.*t20.*t37)-nl3.*(t52-mu.*q21.*t16+kappa.*t3.*t37+mu.*t3.*t20.*t37)-nl1.*(-t39+mu.*q23.*t16+kappa.*t14.*t37+mu.*t14.*t20.*t37);zero;zero;-nl2.*t42+nl1.*t45+nl3.*(mu+kappa.*t46+mu.*t20.*t46);zero;zero;-nl3.*(kappa.*t3.*t49+mu.*t3.*t20.*t49)+nl2.*(-t47+mu.*q31.*t16+kappa.*t9.*t49+mu.*t9.*t20.*t49)-nl1.*(-t51+mu.*q32.*t16+kappa.*t14.*t49+mu.*t14.*t20.*t49);zero;zero;nl3.*(kappa.*t3.*t54+mu.*t3.*t20.*t54)-nl2.*(-t52+mu.*q21.*t16+kappa.*t9.*t54+mu.*t9.*t20.*t54)+nl1.*(-t56+mu.*q22.*t16+kappa.*t14.*t54+mu.*t14.*t20.*t54);zero;zero];
        end
        if nargout > 2
            fb_uh = [-t57;zero;zero;zero;-one;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 8
        kappa = param{2};
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
        t2 = q11.*q32;
        t17 = q12.*q31;
        t3 = t2-t17;
        t4 = q11.*q23.*q32;
        t5 = q12.*q21.*q33;
        t6 = q13.*q22.*q31;
        t10 = q11.*q22.*q33;
        t11 = q12.*q23.*q31;
        t12 = q13.*q21.*q32;
        t7 = t4+t5+t6-t10-t11-t12+1.0;
        t8 = q11.*q33;
        t20 = q13.*q31;
        t9 = t8-t20;
        t13 = q12.*q33;
        t23 = q13.*q32;
        t14 = t13-t23;
        t15 = t4+t5+t6-t10-t11-t12;
        t16 = 1.0./t15;
        fb = [one.*(-uh1+uinf1+x1);-uinf2+tau.*(u2-uh2)+nl3.*(mu.*q23-(mu.*t3)./(t4+t5+t6-q11.*q22.*q33-q12.*q23.*q31-q13.*q21.*q32)+kappa.*t3.*t7)+nl2.*(mu.*q22-kappa.*t7.*t9+mu.*t9.*t16)+nl1.*(mu.*q21+kappa.*t7.*t14-mu.*t14.*t16);-uh3+uinf3+x3];
        if nargout > 1
            t18 = q22.*q33;
            t21 = q23.*q32;
            t19 = t18-t21;
            t22 = 1.0./t15.^2;
            t24 = 1.0./(t4+t5+t6-t10-t11-t12).^2;
            t25 = t14.^2;
            t26 = q12.*q23;
            t28 = q13.*q22;
            t27 = t26-t28;
            t29 = q21.*q33;
            t32 = q23.*q31;
            t30 = t29-t32;
            t31 = mu.*q33.*t16;
            t33 = kappa.*t9.*t14;
            t34 = mu.*t9.*t14.*t24;
            t35 = t33+t34;
            t36 = t9.^2;
            t37 = q11.*q23;
            t39 = q13.*q21;
            t38 = t37-t39;
            t40 = kappa.*q13.*t7;
            t41 = mu.*q31.*t16;
            t42 = q21.*q32;
            t45 = q22.*q31;
            t43 = t42-t45;
            t44 = mu.*q32.*t16;
            t46 = kappa.*t3.*t9;
            t47 = mu.*t3.*t9.*t24;
            t48 = t46+t47;
            t49 = kappa.*t3.*t14;
            t50 = mu.*t3.*t14.*t24;
            t51 = t49+t50;
            t52 = t3.^2;
            t53 = kappa.*q11.*t7;
            t54 = q11.*q22;
            t56 = q12.*q21;
            t55 = t54-t56;
            t57 = kappa.*q12.*t7;
            t58 = one.*tau;
            fb_udg = [zero;zero;zero;zero;t58;zero;zero;zero;zero;zero;-nl1.*(kappa.*t14.*t19+mu.*t14.*t19.*t22)+nl2.*(t31-kappa.*q33.*t7+kappa.*t9.*t19+mu.*t9.*t19.*t22)-nl3.*(t44-kappa.*q32.*t7+kappa.*t3.*t19+mu.*t3.*t19.*t22);zero;zero;-nl2.*t35+nl3.*t51+nl1.*(mu+kappa.*t25+mu.*t24.*t25);zero;zero;-nl1.*(kappa.*t14.*t27+mu.*t14.*t24.*t27)+nl2.*(t40-mu.*q13.*t16+kappa.*t9.*t27+mu.*t9.*t24.*t27)-nl3.*(t57-mu.*q12.*t16+kappa.*t3.*t27+mu.*t3.*t24.*t27);zero;zero;-nl2.*(kappa.*t9.*t30+mu.*t9.*t24.*t30)+nl3.*(t41-kappa.*q31.*t7+kappa.*t3.*t30+mu.*t3.*t24.*t30)+nl1.*(-t31+kappa.*q33.*t7+kappa.*t14.*t30+mu.*t14.*t24.*t30);zero;zero;-nl1.*t35-nl3.*t48+nl2.*(mu+kappa.*t36+mu.*t24.*t36);zero;zero;-nl2.*(kappa.*t9.*t38+mu.*t9.*t24.*t38)+nl3.*(t53-mu.*q11.*t16+kappa.*t3.*t38+mu.*t3.*t24.*t38)+nl1.*(-t40+mu.*q13.*t16+kappa.*t14.*t38+mu.*t14.*t24.*t38);zero;zero;-nl3.*(kappa.*t3.*t43+mu.*t3.*t24.*t43)+nl2.*(-t41+kappa.*q31.*t7+kappa.*t9.*t43+mu.*t9.*t24.*t43)-nl1.*(-t44+kappa.*q32.*t7+kappa.*t14.*t43+mu.*t14.*t24.*t43);zero;zero;-nl2.*t48+nl1.*t51+nl3.*(mu+kappa.*t52+mu.*t24.*t52);zero;zero;-nl3.*(kappa.*t3.*t55+mu.*t3.*t24.*t55)+nl2.*(-t53+mu.*q11.*t16+kappa.*t9.*t55+mu.*t9.*t24.*t55)-nl1.*(-t57+mu.*q12.*t16+kappa.*t14.*t55+mu.*t14.*t24.*t55);zero];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;zero;-t58;zero;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    otherwise
         error('unknown boundary type');
end
