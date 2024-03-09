
fileID = fopen('./run022224_dig/time0084','r');
UDG=fread(fileID,'double');
UDG=reshape(UDG,[6,9,16139]);

Er = UDG(:,6,:);
Ez = UDG(:,9,:);
normE = sqrt(Er.^2+Ez.^2);
figure(1); clf; scaplot(mesh,normE,[],0,0); axis equal; axis tight; axis on; title('|E|');

figure(2); clf; scaplot(mesh,UDG(:,1,:),[],0,0); axis equal; axis tight; axis on; title('ne');

% a=UDG(:,1,:);
% b= a(:);
% c=find(a<0);
% b(c) = .001;
% a = reshape(b,size(UDG(:,1,:)));
% 
% figure(2); clf; scaplot(mesh,log(a),[],0,0); axis equal; axis tight; axis on; title('log(ne)');
