max_loc_vec_2 = zeros(691,1);
max_vec_2 = zeros(691,1);
mesh = mkmesh_streamer_gmsh(2, "streamer_16k-3.msh");

for i=1:400
    disp(i)
    fname = sprintf('run_9-26-23/time%d.mat', i+500);
    load(fname, "UDG");
    [~, Er] = plot_bdry(mesh, UDG(:,:,:), 4, 6);
    [unique_pts, Ez] = plot_bdry(mesh, UDG(:,:,:), 4, 9);

    normE2d = sqrt(Er.^2+Ez.^2)*3e6;
    [M,I] = max(normE2d);
    max_vec_2(i) = M;
    max_loc_vec_2(i) = unique_pts(I,2)*1e-4;
    clear("UDG")
end

mesh = mkmesh_streamer_gmsh(2, "streamer_16k-4.msh");
for i=401:691
    disp(i)
    fname = sprintf('run_9-26-23/time%d.mat', i+500);
    load(fname, "UDG");
    [~, Er] = plot_bdry(mesh, UDG(:,:,:), 4, 6);
    [unique_pts, Ez] = plot_bdry(mesh, UDG(:,:,:), 4, 9);

    normE2d = sqrt(Er.^2+Ez.^2)*3e6;
    [M,I] = max(normE2d);
    max_vec_2(i) = M;
    max_loc_vec_2(i) = unique_pts(I,2)*1e-4;
    clear("UDG")
end
