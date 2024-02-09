itime=170;    % 0-indexed

mesh = mkmesh_streamer_gmsh(porder, "streamer_16k_fixed.msh");
npv=6;
nc=9;
nproc=4;
npf=3;
nch=3;
hybrid='hdg';

[UDG,UH] = getsolfrombinaryfile('./run020824/time',nproc,npv,nc,npf,nch,hybrid, itime);

% fname = sprintf("./run011224/time%04d.bin", itime);
% UDG=readbin(fname);  % udg
% UDG = reshape(UDG, [npv, nc, mesh.ne]);

disp('ne max')
disp(exp(max(max(UDG(:,1,:)))))
figure()
scaplot(mesh,UDG(:,1,:),[],0,0); axis equal; axis tight; colormap jet; title('log(ne)');

disp('E max');
Er = UDG(:,6,:);
Ez = UDG(:,9,:);
normE = sqrt(Er.^2+Ez.^2)*3e6;
max(max(normE))
figure()
scaplot(mesh,normE,[],0,0); axis equal; axis tight; colormap jet; title('|E|');
