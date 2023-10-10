function [f,f_udg] = flux2d(pg,udg,param,time)

cpar = 20.0;
[ng,nc] = size(udg);
nch  = 4;

vis  = param{7};
rmin = param{8};
r = udg(:,1);
ind1 = find(cpar*r/rmin>30);
ind2 = setdiff(1:length(r),ind1);

f = zeros(ng,nch,2);
f_udg = zeros(ng,nch,2,nc);
[f(ind1,:,:),f_udg(ind1,:,:,:)] = flux2d1(pg(ind1,:),udg(ind1,:),param,time);

if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    error('here');
end

[f(ind2,:,:),f_udg(ind2,:,:,:)] = flux2d2(pg(ind2,:),udg(ind2,:),param,time);

if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    error('here');
end

[fav,fav_udg] = avflux2d(pg,udg,param,time);
f = f + vis*bsxfun(@times,pg(:,4),fav);
f_udg = f_udg + vis*bsxfun(@times,pg(:,4),fav_udg);





