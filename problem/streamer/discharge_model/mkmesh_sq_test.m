function mesh = mkmesh_sq_test(porder, fname)

    [p,t] = gmshcall(fname, 2, 0);

    % Normalization beacuse gmsh import is weird for some reason
    xmax = max(p(:,1));
    
    % Below is nondimensional
    bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-6)',...
               'all(p(:,1)>max(p0(:,1))-1e-6)',...
               'all(p(:,2)>max(p0(:,2))-1e-6)',...
               'all(p(:,1)<min(p0(:,1))+1e-6)'};
    
    % Boundaries
    % 1 Bottom electrode
    % 2 Right farfield
    % 3 Top electrode
    % 4 Symmetry
             
    elemtype = 0;
    nodetype = 1;
    mesh = mkmesh(p',t', porder, bndexpr, elemtype, nodetype);