itime=250;

Er = UDG_history(:,6,:,itime);
Ez = UDG_history(:,9,:,itime);
normE = sqrt(Er.^2+Ez.^2);
figure(1); clf; scaplot(mesh,normE*3e6,[],2,0); axis equal; axis tight; axis on;