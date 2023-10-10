

nu=1e-3;
ds = loginc(linspace(0.1,1,500)*4,1);
tetaw = linspace(t1,t2,20); 
x1 = cos(tetaw)*rw+xw;
y1 = sin(tetaw)*rw+yw;

app.arg = {1,nu,eps0,kappa,xw,yw,t1,t2,Wind,tau};
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
VDG(:,1,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,5,:)+Wind*mesh.dgnodes(:,4,:);
[x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x1, y1, ds, 1);

if ~isempty(out)
    Wind2 = Wind;
end
while ~isempty(out)
    Wind = Wind/1.1;    
    app.arg = {1,nu,eps0,kappa,xw,yw,t1,t2,Wind,tau};
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
    VDG(:,1,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,3,:);
    VDG(:,2,:) = UDG(:,5,:)+Wind*mesh.dgnodes(:,4,:);
    [x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x1, y1, ds, 1);
    Wind1 = Wind;
end

if isempty(out)
    Wind1 = Wind;
end
while isempty(out)
    Wind = Wind*1.1;    
    app.arg = {1,nu,eps0,kappa,xw,yw,t1,t2,Wind,tau};
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
    VDG(:,1,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,3,:);
    VDG(:,2,:) = UDG(:,5,:)+Wind*mesh.dgnodes(:,4,:);
    [x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x1, y1, ds, 1);
    Wind2 = Wind;
end

while (1)
    if (Wind2-Wind1)<1e-3 && isempty(out)
        break;
    end    
    Wind = (Wind1+Wind2)/2;    
    app.arg = {1,nu,eps0,kappa,xw,yw,t1,t2,Wind,tau};
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
    VDG(:,1,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,3,:);
    VDG(:,2,:) = UDG(:,5,:)+Wind*mesh.dgnodes(:,4,:);
    [x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x1, y1, ds, 1);    
    if ~isempty(out)
        Wind2 = Wind;
    else
        Wind1 = Wind;
    end
    [Wind1 Wind2 Wind]
end






