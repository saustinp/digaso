
Wind = 0.05;
DV = 35;
Va = -15;
Vw = DV+Va;
t1 = 3.4;
t2 = 5.4;

ds = loginc(linspace(0.1,1,500)*4,1);
tetaw = linspace(t1,t2,20); 
x0 = cos(tetaw)*rw+xw;
y0 = sin(tetaw)*rw+yw;

UDG = UDGa*Va + UDGw*Vw;
VDG(:,1,:) = UDG(:,2,:)+Wind*mesh.dgnodes(:,3,:);
VDG(:,2,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,4,:);
[x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref); 
pause(1)

% all streamlines are attached to the airfoil
if isempty(out)
    Va1 = Va;
end
while isempty(out)
    Va = Va + 1; % increase Va        
    Vw = DV+Va;
    UDG = UDGa*Va + UDGw*Vw;
    VDG(:,1,:) = UDG(:,2,:)+Wind*mesh.dgnodes(:,3,:);
    VDG(:,2,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,4,:);    
    [x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref);
    pause(1)
    Va2 = Va;
end

% Some streamlines are NOT attached to the airfoil
if isempty(out)==0
    Va2 = Va;
end
while isempty(out)==0
    Va = Va - 1; % decrease Va
    Vw = DV+Va;
    UDG = UDGa*Va + UDGw*Vw;
    VDG(:,1,:) = UDG(:,2,:)+Wind*mesh.dgnodes(:,3,:);
    VDG(:,2,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,4,:);    
    [x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref);
    pause(1)
    Va1 = Va;
end
% Va1 < Va2, so Va in [Va1 Va2]

[Va1 Va2]
while (1)
    if abs(Va2-Va1)<1e-2 && isempty(out)
        break;
    end    
    Va = (Va1+Va2)/2;    
    Vw = DV+Va;
    UDG = UDGa*Va + UDGw*Vw;
    VDG(:,1,:) = UDG(:,2,:)+Wind*mesh.dgnodes(:,3,:);
    VDG(:,2,:) = UDG(:,3,:)+Wind*mesh.dgnodes(:,4,:);    
    [x,y,s,Ex,Ey,in,out] = fieldlines(mesh, VDG, x0, y0, ds, nref);    
    pause(1)
    if isempty(out) % all streamlines are attached to the airfoil
        % make Va more positive 
        Va1 = Va; 
    else
        Va2 = Va;
    end
    [Va1 Va Va2]
end






