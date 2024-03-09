Er = UDG(:,6,:);
Ez = UDG(:,9,:);
normE = sqrt(Er.^2+Ez.^2);
figure(1); clf; scaplot(mesh,normE,[],2,0); axis equal; axis tight; axis on;
figure(2); clf; scaplot(mesh,UDG(:,1,:),[],2,0); colormap jet; axis equal; axis tight; axis on;