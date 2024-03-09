function mesh = mkmesh_dual_electrode(porder, fname)

    [p,t] = gmshcall(fname, 2, 0);

    % Normalization beacuse gmsh import is weird for some reason
    xmax = max(p(:,1));
    % p = p/xmax * 100;
    
    % Below is nondimensional
    bndexpr = {'all(p(:,1)<min(p0(:,1))+1e-4)',...
               'all(p(:,1)>max(p0(:,1))-1e-4) | all(p(:,2)>max(p0(:,2))-1e-4) | all(p(:,2)<min(p0(:,2))+1e-4)',...
                'all(p(:,2)>-20)',...
                'all(p(:,2)<-20)'};
    
    % Boundaries
    % 1 Symmetry
    % 2 Right farfield
    % 3 Top electrode
    % 4 Bottom electrode
             
    elemtype = 0;
    nodetype = 1;
    mesh = mkmesh(p',t', porder, bndexpr, elemtype, nodetype);