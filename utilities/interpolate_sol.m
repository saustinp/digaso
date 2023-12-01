function UDG = interpolate_sol(mesh1, mesh2, UDG1)

npe=size(mesh2.dgnodes,1);
nc=size(UDG1,2);
[elist, xi] = locatexinmesh(mesh1, mesh2.dgnodes, [], 1e-4);
UDG   = evalfield(mesh1, UDG1, elist, xi);
UDG = permute(reshape(UDG, [npe mesh2.ne nc]),[1 3 2]); 

end