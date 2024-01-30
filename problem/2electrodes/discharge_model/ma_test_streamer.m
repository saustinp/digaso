load('/Users/saustin/Documents/HDG/problem/2d_streamer/discharge_model/run_9-26-23/time1000.mat')
mesh1 = mkmesh_streamer_gmsh(2, "streamer_16k-4.msh");  % Original mesh that the solution was computed on

porder=2;
mesh2 = mkmesh_streamer(50,50,2,0,[0 125 0 125],0,1);   % Coarser mesh to interpolate the solution onto for testing MA
master = mkmaster(mesh2,2*porder);
[master,mesh2] = preprocess(master,mesh2,hybrid);
mesh2 = mkcgmesh(mesh2);
meshplot(mesh2);

UDG_coarse = interpolate_sol(mesh1, mesh2, UDG);

% Plot the normE
Er = UDG_coarse(:,6,:);
Ez = UDG_coarse(:,9,:);
normE = sqrt(Er.^2+Ez.^2)*3e6;  % This will be used as the sensor field

normE = normE/max(max(normE))*10;   % Rescaling to [0,10]

scaplot(mesh2,normE,[],0,1); axis equal; axis tight; colormap jet;

