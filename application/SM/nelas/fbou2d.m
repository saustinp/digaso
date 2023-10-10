function [fb,fb_udg,fb_uh] = fbou2d(ib,uinf,nl,pg,udg,uh,param,time)
%FBOU2D.M
%    [FB,FB_UDG,FB_UH] = FBOU2D(IB,UINF,NL,PG,UDG,UH,PARAM,TIME)
[ng,nc] = size(udg);
nch = 2;
nd = 2;
switch (ib)
    case 1
        one = ones(ng,1);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        x1 = pg(:,1);
        x2 = pg(:,2);
        zero = zeros(ng,1);
        fb = [one.*(-uh1+uinf1+x1);-uh2+uinf2+x2];
        if nargout > 1
            fb_udg = [zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 2
        [fb,fb_udg,fb_uh] = fhat(nl,pg,udg,uh,param,time);
        fb = fb - uinf;
    case 3
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        one = ones(ng,1);
        pres = udg(:,3);
        q11 = udg(:,4);
        q12 = udg(:,6);
        q21 = udg(:,5);
        q22 = udg(:,7);
        tau = param{3};
        u2 = udg(:,2);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        x1 = pg(:,1);
        zero = zeros(ng,1);
        t2032 = q11.*q22;
        t2035 = q12.*q21;
        t2033 = t2032-t2035;
        t2034 = 1.0./t2033;
        fb = [one.*(-uh1+uinf1+x1);-uinf2+nl1.*(mu.*q21-pres.*q12+mu.*q12.*t2034)+nl2.*(mu.*q22+pres.*q11-mu.*q11.*t2034)+tau.*(u2-uh2)];
        if nargout > 1
            t2036 = 1.0./t2033.^2;
            t2037 = one.*tau;
            fb_udg = [zero;zero;zero;t2037;zero;-nl1.*q12+nl2.*q11;zero;nl2.*(pres-mu.*t2034+mu.*q11.*q22.*t2036)-mu.*nl1.*q12.*q22.*t2036;zero;nl1.*(mu+mu.*q12.^2.*t2036)-mu.*nl2.*q11.*q12.*t2036;zero;nl1.*(-pres+mu.*t2034+mu.*q12.*q21.*t2036)-mu.*nl2.*q11.*q21.*t2036;zero;nl2.*(mu+mu.*q11.^2.*t2036)-mu.*nl1.*q11.*q12.*t2036];
        end
        if nargout > 2
            fb_uh = [-one;zero;zero;-t2037];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    case 4
        mu = param{1};
        nl1 = nl(:,1);
        nl2 = nl(:,2);
        one = ones(ng,1);
        pres = udg(:,3);
        q11 = udg(:,4);
        q12 = udg(:,6);
        q21 = udg(:,5);
        q22 = udg(:,7);
        tau = param{3};
        u1 = udg(:,1);
        uh1 = uh(:,1);
        uh2 = uh(:,2);
        uinf1 = uinf(:,1);
        uinf2 = uinf(:,2);
        x2 = pg(:,2);
        zero = zeros(ng,1);
        t2039 = q11.*q22;
        t2042 = q12.*q21;
        t2040 = -t2042+t2039;
        t2041 = 1.0./t2040;
        fb = [-uinf1+nl1.*(mu.*q11+pres.*q22-mu.*q22.*t2041)+nl2.*(mu.*q12-pres.*q21+mu.*q21.*t2041)+tau.*(u1-uh1);-uh2+uinf2+x2];
        if nargout > 1
            t2043 = 1.0./t2040.^2;
            t2044 = mu.*t2041;
            t2045 = one.*tau;
            fb_udg = [t2045;zero;zero;zero;nl1.*q22-nl2.*q21;zero;nl1.*(mu+mu.*q22.^2.*t2043)-mu.*nl2.*q21.*q22.*t2043;zero;nl2.*(-pres+t2044+mu.*q12.*q21.*t2043)-mu.*nl1.*q12.*q22.*t2043;zero;nl2.*(mu+mu.*q21.^2.*t2043)-mu.*nl1.*q21.*q22.*t2043;zero;nl1.*(pres-t2044+mu.*q11.*q22.*t2043)-mu.*nl2.*q11.*q21.*t2043;zero];
        end
        if nargout > 2
            fb_uh = [-t2045;zero;zero;-one];
        end
        fb = reshape(fb,ng,nch);
        fb_udg = reshape(fb_udg,ng,nch,nc);
        fb_uh = reshape(fb_uh,ng,nch,nch);
    otherwise
         error('unknown boundary type');
end
