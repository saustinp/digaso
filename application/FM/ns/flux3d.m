function [f,f_udg] = flux3d(pg,udg,param,time)
%FLUX3D
%    [F,F_UDG] = FLUX3D(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    06-Nov-2012 11:00:33
[ng,nc] = size(udg);
nch = 5;
nd = 3;
one = ones(ng,1);
param1 = param{1};
param3 = param{3};
param4 = param{4};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u4 = udg(:,4);
u5 = udg(:,5);
u6 = udg(:,6);
u7 = udg(:,7);
u8 = udg(:,8);
u9 = udg(:,9);
u10 = udg(:,10);
u11 = udg(:,11);
u12 = udg(:,12);
u13 = udg(:,13);
u14 = udg(:,14);
u15 = udg(:,15);
u16 = udg(:,16);
u17 = udg(:,17);
u18 = udg(:,18);
u19 = udg(:,19);
u20 = udg(:,20);
zero = zeros(ng,1);
t274 = 1.0./u1;
t275 = u2.^2;
t276 = 1.0./u1.^2;
t277 = 1.0./param3;
t278 = t275.*t276.*(1.0./2.0);
t279 = u3.^2;
t280 = t276.*t279.*(1.0./2.0);
t281 = u4.^2;
t282 = t276.*t281.*(1.0./2.0);
t283 = t278+t280+t282;
t305 = t283.*u1;
t284 = -t305+u5;
t285 = param1-1.0;
t303 = t274.*u3.*u6;
t286 = -t303+u8;
t287 = t274.*t286;
t306 = t274.*u2.*u11;
t288 = -t306+u12;
t289 = t274.*t288;
t290 = t287+t289;
t304 = t274.*u4.*u6;
t291 = -t304+u9;
t292 = t274.*t291;
t326 = t274.*u2.*u16;
t293 = -t326+u17;
t294 = t274.*t293;
t295 = t292+t294;
t302 = t274.*u2.*u6;
t296 = -t302+u7;
t310 = t274.*u3.*u11;
t297 = -t310+u13;
t298 = t274.*t297;
t311 = t274.*u4.*u16;
t299 = -t311+u19;
t300 = t274.*t299;
t341 = t274.*t296.*2.0;
t301 = t298+t300-t341;
t307 = t277.*t290;
t308 = t274.*u2.*u3;
t309 = t307+t308;
t312 = t284.*t285;
t313 = t274.*u5;
t314 = t274.*t284.*t285;
t315 = t313+t314;
t324 = t274.*u4.*u11;
t316 = -t324+u14;
t317 = t274.*t316;
t330 = t274.*u3.*u16;
t318 = -t330+u18;
t319 = t274.*t318;
t320 = t317+t319;
t321 = t274.*t296;
t372 = t274.*t297.*2.0;
t322 = t300+t321-t372;
t323 = 1.0./param4;
t325 = 1.0./t285;
t327 = t277.*t295;
t328 = t274.*u2.*u4;
t329 = t327+t328;
t331 = t277.*t320;
t332 = t274.*u3.*u4;
t333 = t331+t332;
t397 = t274.*t299.*2.0;
t334 = t298+t321-t397;
t335 = 1.0./u1.^3;
t336 = t275.*t335;
t337 = t279.*t335;
t338 = t281.*t335;
t339 = t336+t337+t338;
t359 = t339.*u1;
t340 = t278+t280+t282-t359;
t342 = t276.*t286;
t343 = t276.*t288;
t364 = t335.*u3.*u6;
t365 = t335.*u2.*u11;
t344 = t342+t343-t364-t365;
t345 = t276.*t291;
t346 = t276.*t293;
t389 = t335.*u4.*u6;
t390 = t335.*u2.*u16;
t347 = t345+t346-t389-t390;
t348 = t276.*t296.*2.0;
t349 = t335.*u3.*u11;
t350 = t335.*u4.*u16;
t367 = t276.*t299;
t396 = t276.*t297;
t351 = t348+t349+t350-t367-t396-t335.*u2.*u6.*2.0;
t352 = t276.*t296.*u2;
t353 = t276.*t286.*u3;
t354 = t276.*t291.*u4;
t355 = t352+t353+t354;
t356 = t355.*u1;
t357 = t283.*u6;
t358 = t356+t357-u10;
t360 = 1.0./u1.^4;
t361 = t285.*t358.*u1;
t362 = t284.*t285.*u6;
t363 = t361+t362;
t366 = -t277.*t344-t276.*u2.*u3;
t368 = t276.*u5;
t369 = t276.*t284.*t285;
t370 = t274.*t285.*t340;
t371 = t368+t369+t370;
t373 = t276.*t316;
t374 = t276.*t318;
t392 = t335.*u4.*u11;
t393 = t335.*u3.*u16;
t375 = t373+t374-t392-t393;
t376 = t276.*t297.*2.0;
t377 = t335.*u2.*u6;
t395 = t276.*t296;
t378 = t350-t367+t376+t377-t395-t335.*u3.*u11.*2.0;
t379 = t276.*t288.*u2;
t380 = t276.*t297.*u3;
t381 = t276.*t316.*u4;
t382 = t379+t380+t381;
t383 = t382.*u1;
t384 = t283.*u11;
t385 = t383+t384-u15;
t386 = t285.*t385.*u1;
t387 = t284.*t285.*u11;
t388 = t386+t387;
t391 = -t277.*t347-t276.*u2.*u4;
t394 = -t277.*t375-t276.*u3.*u4;
t398 = t276.*t299.*2.0;
t399 = t349+t377-t395-t396+t398-t335.*u4.*u16.*2.0;
t400 = t276.*t293.*u2;
t401 = t276.*t318.*u3;
t402 = t276.*t299.*u4;
t403 = t400+t401+t402;
t404 = t403.*u1;
t405 = t283.*u16;
t406 = t404+t405-u20;
t407 = t285.*t406.*u1;
t408 = t284.*t285.*u16;
t409 = t407+t408;
f = [u2;t312+t274.*t275-t277.*t301.*(2.0./3.0);t309;t329;t315.*u2+t274.*t277.*t290.*u3+t274.*t277.*t295.*u4-t274.*t277.*t301.*u2.*(2.0./3.0)-param1.*t276.*t277.*t323.*t325.*t363;u3;t309;t312+t274.*t279-t277.*t322.*(2.0./3.0);t333;t315.*u3+t274.*t277.*t290.*u2+t274.*t277.*t320.*u4-t274.*t277.*t322.*u3.*(2.0./3.0)-param1.*t276.*t277.*t323.*t325.*t388;u4;t329;t333;t312+t274.*t281-t277.*t334.*(2.0./3.0);t315.*u4+t274.*t277.*t295.*u2+t274.*t277.*t320.*u3-t274.*t277.*t334.*u4.*(2.0./3.0)-param1.*t276.*t277.*t323.*t325.*t409];
if nargout > 1
    t410 = t274.*u3;
    t426 = t276.*t277.*u11;
    t411 = t410-t426;
    t412 = t274.*u4;
    t419 = t276.*t277.*u16;
    t413 = t412-t419;
    t414 = t276.*t277.*u6.*(2.0./3.0);
    t415 = t414-t274.*t285.*u2;
    t416 = t274.*t277.*t290;
    t417 = t274.*u2;
    t422 = t276.*t277.*u6;
    t418 = t417-t422;
    t420 = t276.*t277.*u11.*(2.0./3.0);
    t421 = t420-t274.*t285.*u3;
    t423 = t274.*t277.*t295;
    t424 = t276.*t277.*u16.*(2.0./3.0);
    t428 = t274.*t285.*u4;
    t425 = t424-t428;
    t427 = t274.*t277.*t320;
    t429 = one.*t285;
    t430 = t274.*t285;
    t431 = t274+t430;
    t432 = t276.*t277.*u2.*(2.0./3.0);
    t433 = t274.*t277;
    t434 = t276.*t277.*u2;
    t435 = t276.*t277.*u4;
    t436 = t285.*t340.*u1;
    t437 = t312+t436;
    t438 = t276.*t277.*u3.*(2.0./3.0);
    t439 = t276.*t277.*u3;
    t440 = t274.*t277.*(4.0./3.0);
    t446 = param1.*t276.*t277.*t323.*u4;
    t441 = t435-t446;
    t442 = param1.*t274.*t277.*t323;
    t443 = t276.*t277.*u4.*(2.0./3.0);
    t444 = t434-param1.*t276.*t277.*t323.*u2;
    t445 = t439-param1.*t276.*t277.*t323.*u3;
    f_udg = [zero;-t275.*t276-t285.*t340-t277.*t351.*(2.0./3.0);t366;t391;-t371.*u2-t276.*t277.*t290.*u3-t276.*t277.*t295.*u4+t276.*t277.*t301.*u2.*(2.0./3.0)-t274.*t277.*t344.*u3-t274.*t277.*t347.*u4-t274.*t277.*t351.*u2.*(2.0./3.0)+param1.*t277.*t323.*t325.*t335.*t363.*2.0-param1.*t276.*t277.*t323.*t325.*(t285.*t358-t285.*t340.*u6+t285.*u1.*(t352+t353+t354-t339.*u6-u1.*(t286.*t335.*u3.*2.0+t291.*t335.*u4.*2.0+t296.*t335.*u2.*2.0-t275.*t360.*u6-t279.*t360.*u6-t281.*t360.*u6)));zero;t366;-t276.*t279-t285.*t340-t277.*t378.*(2.0./3.0);t394;-t371.*u3-t276.*t277.*t290.*u2-t276.*t277.*t320.*u4+t276.*t277.*t322.*u3.*(2.0./3.0)-t274.*t277.*t344.*u2-t274.*t277.*t375.*u4-t274.*t277.*t378.*u3.*(2.0./3.0)+param1.*t277.*t323.*t325.*t335.*t388.*2.0-param1.*t276.*t277.*t323.*t325.*(t285.*t385-t285.*t340.*u11+t285.*u1.*(t379+t380+t381-t339.*u11-u1.*(t288.*t335.*u2.*2.0+t297.*t335.*u3.*2.0-t275.*t360.*u11-t279.*t360.*u11-t281.*t360.*u11+t316.*t335.*u4.*2.0)));zero;t391;t394;-t276.*t281-t285.*t340-t277.*t399.*(2.0./3.0);-t371.*u4-t276.*t277.*t295.*u2-t276.*t277.*t320.*u3+t276.*t277.*t334.*u4.*(2.0./3.0)-t274.*t277.*t347.*u2-t274.*t277.*t375.*u3-t274.*t277.*t399.*u4.*(2.0./3.0)+param1.*t277.*t323.*t325.*t335.*t409.*2.0-param1.*t276.*t277.*t323.*t325.*(t285.*t406-t285.*t340.*u16+t285.*u1.*(t400+t401+t402-t339.*u16-u1.*(t293.*t335.*u2.*2.0+t299.*t335.*u4.*2.0-t275.*t360.*u16-t279.*t360.*u16+t318.*t335.*u3.*2.0-t281.*t360.*u16)));one;t274.*u2.*2.0-t276.*t277.*u6.*(4.0./3.0)-t274.*t285.*u2;t411;t413;t313+t314-t275.*t276.*t285-t274.*t277.*t301.*(2.0./3.0)-t277.*t335.*u2.*u6.*(4.0./3.0)-t277.*t335.*u3.*u11-t277.*t335.*u4.*u16+param1.*t276.*t277.*t323.*t325.*(t285.*u1.*(u1.*(t377-t395)-t276.*u2.*u6)+t274.*t285.*u2.*u6);zero;t411;t415;zero;t416-t276.*t285.*u2.*u3+t277.*t335.*u3.*u6.*(2.0./3.0)-t277.*t335.*u2.*u11-param1.*t276.*t277.*t323.*t325.*(t285.*u1.*(u1.*(t343-t365)+t276.*u2.*u11)-t274.*t285.*u2.*u11);zero;t413;zero;t415;t423-t276.*t285.*u2.*u4+t277.*t335.*u4.*u6.*(2.0./3.0)-t277.*t335.*u2.*u16-param1.*t276.*t277.*t323.*t325.*(t285.*u1.*(u1.*(t346-t390)+t276.*u2.*u16)-t274.*t285.*u2.*u16);zero;t421;t418;zero;t416-t276.*t285.*u2.*u3-t277.*t335.*u3.*u6+t277.*t335.*u2.*u11.*(2.0./3.0)-param1.*t276.*t277.*t323.*t325.*(t285.*u1.*(u1.*(t342-t364)+t276.*u3.*u6)-t274.*t285.*u3.*u6);one;t418;t274.*u3.*2.0-t274.*t285.*u3-t276.*t277.*u11.*(4.0./3.0);t413;t313+t314-t276.*t279.*t285-t274.*t277.*t322.*(2.0./3.0)-t277.*t335.*u2.*u6-t277.*t335.*u3.*u11.*(4.0./3.0)-t277.*t335.*u4.*u16+param1.*t276.*t277.*t323.*t325.*(t285.*u1.*(u1.*(t349-t396)-t276.*u3.*u11)+t274.*t285.*u3.*u11);zero;zero;t413;t421;t427-t276.*t285.*u3.*u4+t277.*t335.*u4.*u11.*(2.0./3.0)-t277.*t335.*u3.*u16-param1.*t276.*t277.*t323.*t325.*(t285.*u1.*(u1.*(t374-t393)+t276.*u3.*u16)-t274.*t285.*u3.*u16);zero;t425;zero;t418;t423-t276.*t285.*u2.*u4-t277.*t335.*u4.*u6+t277.*t335.*u2.*u16.*(2.0./3.0)-param1.*t276.*t277.*t323.*t325.*(t285.*u1.*(u1.*(t345-t389)+t276.*u4.*u6)-t274.*t285.*u4.*u6);zero;zero;t425;t411;t427-t276.*t285.*u3.*u4-t277.*t335.*u4.*u11+t277.*t335.*u3.*u16.*(2.0./3.0)-param1.*t276.*t277.*t323.*t325.*(t285.*u1.*(u1.*(t373-t392)+t276.*u4.*u11)-t274.*t285.*u4.*u11);one;t418;t411;-t428+t274.*u4.*2.0-t276.*t277.*u16.*(4.0./3.0);t313+t314-t276.*t281.*t285-t274.*t277.*t334.*(2.0./3.0)-t277.*t335.*u2.*u6-t277.*t335.*u3.*u11-t277.*t335.*u4.*u16.*(4.0./3.0)+param1.*t276.*t277.*t323.*t325.*(t285.*u1.*(u1.*(t350-t367)-t276.*u4.*u16)+t274.*t285.*u4.*u16);zero;t429;zero;zero;t431.*u2-param1.*t276.*t277.*t323.*u6;zero;zero;t429;zero;t431.*u3-param1.*t276.*t277.*t323.*u11;zero;zero;zero;t429;t431.*u4-param1.*t276.*t277.*t323.*u16;zero;t276.*t277.*u2.*(-4.0./3.0);-t276.*t277.*u3;-t276.*t277.*u4;t275.*t277.*t335.*(-4.0./3.0)-t277.*t279.*t335-t277.*t281.*t335-param1.*t276.*t277.*t323.*t325.*t437;zero;-t276.*t277.*u3;t432;zero;t277.*t335.*u2.*u3.*(-1.0./3.0);zero;-t276.*t277.*u4;zero;t432;t277.*t335.*u2.*u4.*(-1.0./3.0);zero;t440;zero;zero;t276.*t277.*u2.*(4.0./3.0)-param1.*t276.*t277.*t323.*u2;zero;zero;t274.*t277.*(-2.0./3.0);zero;t276.*t277.*u3.*(-2.0./3.0);zero;zero;zero;t274.*t277.*(-2.0./3.0);t276.*t277.*u4.*(-2.0./3.0);zero;zero;t433;zero;t445;zero;t433;zero;zero;t434;zero;zero;zero;zero;zero;zero;zero;zero;t433;t441;zero;zero;zero;zero;zero;zero;t433;zero;zero;t434;zero;zero;zero;zero;t442;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;t438;-t434;zero;t277.*t335.*u2.*u3.*(-1.0./3.0);zero;-t434;t276.*t277.*u3.*(-4.0./3.0);-t435;-t275.*t277.*t335-t277.*t279.*t335.*(4.0./3.0)-t277.*t281.*t335-param1.*t276.*t277.*t323.*t325.*t437;zero;zero;-t435;t438;t277.*t335.*u3.*u4.*(-1.0./3.0);zero;zero;t433;zero;t439;zero;t433;zero;zero;t444;zero;zero;zero;zero;zero;zero;t274.*t277.*(-2.0./3.0);zero;zero;-t432;zero;zero;t440;zero;t276.*t277.*u3.*(4.0./3.0)-param1.*t276.*t277.*t323.*u3;zero;zero;zero;t274.*t277.*(-2.0./3.0);t276.*t277.*u4.*(-2.0./3.0);zero;zero;zero;zero;zero;zero;zero;zero;t433;t441;zero;zero;t433;zero;t439;zero;zero;zero;zero;zero;zero;zero;zero;zero;t442;zero;zero;zero;zero;zero;zero;t443;zero;-t434;t277.*t335.*u2.*u4.*(-1.0./3.0);zero;zero;t443;-t439;t277.*t335.*u3.*u4.*(-1.0./3.0);zero;-t434;-t439;t276.*t277.*u4.*(-4.0./3.0);-t275.*t277.*t335-t277.*t279.*t335-t277.*t281.*t335.*(4.0./3.0)-param1.*t276.*t277.*t323.*t325.*t437;zero;zero;zero;t433;t435;zero;zero;zero;zero;zero;zero;t433;zero;zero;t444;zero;zero;zero;zero;zero;zero;zero;zero;t433;t435;zero;zero;t433;zero;t445;zero;t274.*t277.*(-2.0./3.0);zero;zero;-t432;zero;zero;t274.*t277.*(-2.0./3.0);zero;-t438;zero;zero;zero;t440;-t446+t276.*t277.*u4.*(4.0./3.0);zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;t442];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);