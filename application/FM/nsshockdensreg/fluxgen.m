syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20
syms x1 x2 x3 
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10
syms zero one 

param = [param1 param2 param3 param4 param5 param6 param7 param8 param9 param10];
if nd==2    
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12];        
    pg = [x1 x2];       
else
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20];        
    pg = [x1 x2 x3];      
end

gam  = param(1);
gam1 = param(1) - 1.0;          
Re   = param(3);
Pr   = param(4);
Minf = param(5);
tau  = param(6);
vis  = param(7);
rmin = param(8);
Re1  = 1/Re;
M2   = Minf^2;
c23  = 2/3;
fc = 1/(gam1*M2*Re*Pr);

hpar = 1.0;
apar = 10.0;
bpar = 0.5;
cpar = 20.0;

ncu = nd+2;
nc = length(udg);
nch = ncu;

if nd==2                                               
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rE   = udg(4);
    rx   = udg(5);
    rux  = udg(6);
    rvx  = udg(7);
    rEx  = udg(8);
    ry   = udg(9);
    ruy  = udg(10);
    rvy  = udg(11);
    rEy  = udg(12);

    r1   = 1/r;
    u    = ru*r1;
    v    = rv*r1;
    E    = rE*r1;
    q    = 0.5*(u*u+v*v);
    p    = gam1*(rE-r*q);
    h    = E+p*r1;    
    
    fi   = [ru, ru*u+p, rv*u,   ru*h, ...
            rv, ru*v,   rv*v+p, rv*h];            
        
    ux  = (rux - rx*u)*r1;
    vx  = (rvx - rx*v)*r1;
    qx  = u*ux + v*vx;
    px  = gam1*(rEx - rx*q - r*qx);
    Tx  = gam*M2*(px*r - p*rx)*r1^2;

    uy  = (ruy - ry*u)*r1;
    vy  = (rvy - ry*v)*r1;
    qy  = u*uy + v*vy;
    py  = gam1*(rEy - ry*q - r*qy);
    Ty  = gam*M2*(py*r - p*ry)*r1^2;

    txx = Re1*c23*(2*ux - vy);
    txy = Re1*(uy + vx);
    tyy = Re1*c23*(2*vy - ux);

    fv = [0, txx, txy, u*txx + v*txy + fc*Tx,...
          0, txy, tyy, u*txy + v*tyy + fc*Ty];        
      
    f1 = fi+fv;  
    
    rt   = (1/cpar)*log(cpar+exp(cpar*(r/rmin)))*rmin;
    r1   = 1/rt;
    u    = ru*r1;
    v    = rv*r1;
    E    = rE*r1;
    q    = 0.5*(u*u+v*v);
    p    = gam1*(rE-rt*q);
    h    = E+p*r1;    
    
    fi   = [ru, ru*u+p, rv*u,   ru*h, ...
            rv, ru*v,   rv*v+p, rv*h];            
        
    rtx = rx*(exp((cpar*r)/rmin)/(cpar+exp((cpar*r)/rmin)));
    ux  = (rux - rtx*u)*r1;
    vx  = (rvx - rtx*v)*r1;
    qx  = u*ux + v*vx;
    px  = gam1*(rEx - rtx*q - rt*qx);
    Tx  = gam*M2*(px*rt - p*rtx)*r1^2;

    rty = ry*(exp((cpar*r)/rmin)/(cpar+exp((cpar*r)/rmin)));
    uy  = (ruy - rty*u)*r1;
    vy  = (rvy - rty*v)*r1;
    qy  = u*uy + v*vy;
    py  = gam1*(rEy - rty*q - rt*qy);
    Ty  = gam*M2*(py*rt - p*rty)*r1^2;

    txx = Re1*c23*(2*ux - vy);
    txy = Re1*(uy + vx);
    tyy = Re1*c23*(2*vy - ux);

    fv = [0, txx, txy, u*txx + v*txy + fc*Tx,...
          0, txy, tyy, u*txy + v*tyy + fc*Ty];            
      
    f2 = fi+fv;        
else
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rw   = udg(4);
    rE   = udg(5);

    rx   = udg(6);
    rux  = udg(7);
    rvx  = udg(8);
    rwx  = udg(9);
    rEx  = udg(10);

    ry   = udg(11);
    ruy  = udg(12);
    rvy  = udg(13);
    rwy  = udg(14);
    rEy  = udg(15);    
end

for n=1:2
    if n==1
        f = f1;
        filename1 = ['flux' num2str(nd) 'd1' '.m'];
    elseif n==2
        f = f2;
        filename1 = ['flux' num2str(nd) 'd2' '.m'];
    end
    
    %%% compute Jacobian
    jac_f  = jacobian(f,udg);

    %%% And patch with vector zero or one to have the right sizes
    for ii = 1:size(f,1)
        for jj = 1:size(f,2)
            temp = ismember(symvar(f(ii,jj)),udg);
            if f(ii,jj)==0, f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, f(ii,jj) = (f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

    for ii = 1:size(jac_f,1)
        for jj = 1:size(jac_f,2)
            temp = ismember(symvar(jac_f(ii,jj)),udg);
            if jac_f(ii,jj)==0, jac_f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, jac_f(ii,jj) = (jac_f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end    

    % generate a temporary matlab file
    matlabFunction(f(:),jac_f(:),'file','tmp.m',...
        'vars',{pg,udg,param,time,[zero one]},'outputs', {'f','f_udg'});

    % open the file and modify it
    fid = fopen('tmp.m','r');
    gid = fopen(filename1,'wt');

    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        str = strrep(tline, 'tmp', strrep(filename1,'.m',''));
        str = strrep(str, 'TMP', upper(strrep(filename1,'.m','')));
        str = strrep(str, 'in1', 'pg');
        str = strrep(str, 'IN1', 'PG');        
        str = strrep(str, 'in2', 'udg');
        str = strrep(str, 'IN2', 'UDG');        
        str = strrep(str, 'in3', 'param');
        str = strrep(str, 'IN3', 'PARAM');                
        str = strrep(str, ',in5)', ')');                
        str = strrep(str, ',IN5)', ')');    
        str = strrep(str, 'param(:,1)', 'param{1}'); 
        str = strrep(str, 'param(:,2)', 'param{2}'); 
        str = strrep(str, 'param(:,3)', 'param{3}'); 
        str = strrep(str, 'param(:,4)', 'param{4}'); 
        str = strrep(str, 'param(:,5)', 'param{5}'); 
        str = strrep(str, 'param(:,6)', 'param{6}'); 
        str = strrep(str, 'param(:,7)', 'param{7}'); 
        str = strrep(str, 'param(:,8)', 'param{8}'); 
        str = strrep(str, 'param(:,9)', 'param{9}');     
        str = strrep(str, 'param(:,10)', 'param{10}');     
        str = strrep(str, 'in5(:,1)', 'zeros(ng,1)');
        str = strrep(str, 'in5(:,2)', 'ones(ng,1)');        
        if i==7
            str = '[ng,nc] = size(udg);';
            fprintf(gid, '%s\n', str);                  
            str = ['nch = ' num2str(nch) ';'];
            fprintf(gid, '%s\n', str);                  
            str = ['nd = ' num2str(nd) ';'];
        end
        fprintf(gid, '%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;
        %disp(str)
    end

    str = 'f = reshape(f,ng,nch,nd);';
    fprintf(gid, '%s\n', str);                  
    str = 'f_udg = reshape(f_udg,ng,nch,nd,nc);';
    fprintf(gid, '%s\n', str);                  

    fclose(fid);
    fclose(gid);
    delete('tmp.m');
end

