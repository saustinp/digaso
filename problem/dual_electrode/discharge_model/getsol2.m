% % mesh = mkmesh_dual_electrode(2, "delec40k.msh");
% [UDG,~] = getsolfrombinaryfile('./run022724_d1/soltime', './run022724_d1/streamersol',6,6,9,3,3,1, 10);   
% Er = UDG(:,6,:);
% Ez = UDG(:,9,:);
% normE = sqrt(Er.^2+Ez.^2)*3e6;
% figure(); scaplot(mesh,normE,[],0,0); axis equal; axis tight; colormap jet; title('|E|');

% scaplot(mesh,UDG(:,1,:),[],0,1); axis equal; axis tight; colormap jet; title('|E|');

[UDG,~] = getsolfrombinaryfile('./run022724_d1/soltime', './run022724_d1/streamersol',6,6,9,3,3,1, 300);   
Er = UDG(:,6,:);
Ez = UDG(:,9,:);
normE = sqrt(Er.^2+Ez.^2)*3e6;
figure(); scaplot(mesh,normE,[],0,0); axis equal; axis tight; colormap jet; title('|E|');