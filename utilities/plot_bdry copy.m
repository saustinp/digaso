function [unique_pts, bdry_sol] = plot_bdry(mesh, UDG, bdry_ind, sol_field)

if size(UDG) ~= mesh.ne
    error("number of elements and size of UDG do not match, perhaps you're not using the right mesh?")
end

bdry_face_idx = find(mesh.f(:,4) == -bdry_ind);
elems_on_bdry = mesh.f(bdry_face_idx,3);
ndim = size(mesh.dgnodes,2);
nodesPerFace = size(mesh.perm,1);

dgSol = zeros(nodesPerFace*size(bdry_face_idx,1),1);
dgPts = zeros(nodesPerFace*size(bdry_face_idx,1),ndim);

for i=1:size(elems_on_bdry)
    elem_idx = elems_on_bdry(i);
    locFaceInd = find(mesh.t2f(elem_idx,:)==bdry_face_idx(i));
    locNodes = mesh.perm(:,locFaceInd);
    if any(UDG~=0, 'all')  % Accepts the case in which UDG is not passed in
        dgSol(1+(i-1)*nodesPerFace:i*nodesPerFace) = UDG(locNodes,sol_field,elem_idx);
    end
    dgPts(1+(i-1)*nodesPerFace:i*nodesPerFace, :) = mesh.dgnodes(locNodes,:,elem_idx);
end

% Remove nonunique points and values from the solution vector
snap=1e-8;
A=round(dgPts/snap)*snap;

[unique_pts,unique_ind] = unique(A, 'rows');
bdry_sol = dgSol(unique_ind);
