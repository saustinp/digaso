% porder=2;
% mesh1 = mkmesh_rect(15,200,porder,0,[0 5 0 125],0,1);
% bdry_pts = boundarypoints(mesh1.p,mesh1.f,2);
% pts2 = [bdry_pts; 125 125; 125 0;];
% pts2 = pts2(end:-1:1,:);
% poly2gmsh('rightdomain.geo', pts2, 5);
% disp('Run the command in a terminal: gmsh rightdomain.geo -2 -o rightdomain.msh3 and then rename to remove the 3');
% return;

[gmshp, gmsht] = gmshcall('rightdomain.msh', 2, 0);
gmshp = gmshp';
gmsht = gmsht';

[p,t] = connectmesh(mesh1.p,mesh1.t,gmshp,gmsht,1e-2);

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-6)',...
           'all(p(:,1)>max(p0(:,1))-1e-6)',...
           'all(p(:,2)>max(p0(:,2))-1e-6)',...
           'all(p(:,1)<min(p0(:,1))+1e-6)'};

mesh2 = mkmesh(p,t,mesh.porder,bndexpr,0,1);
figure(1);clf; meshplot(mesh2);