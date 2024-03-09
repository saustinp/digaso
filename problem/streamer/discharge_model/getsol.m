% function[UDG,UH] = getsol(timestep)

time = timestep;
nproc = 4;
filename = "run021424/time";

elempart = []; entpart = [];
nelem = zeros(nproc,1);
nent = zeros(nproc,1);
for i = 1:nproc
    fileID = fopen([filename + 'elempart_np' + num2str(i-1) + '.bin'],'r');
    elemparti = fread(fileID,'int');
    nelem(i) = length(elemparti);
    elempart = [elempart; elemparti];
    fclose(fileID);

    fileID = fopen([filename + 'entpart_np' + num2str(i-1) + '.bin'],'r');
    entparti = fread(fileID,'int');
    nent(i) = length(entparti);
    entpart = [entpart; entparti];
    fclose(fileID);
end
elempart = elempart + 1;
entpart = entpart+1;

n1 = cumsum([0; nelem]); 
n2 = cumsum([0; nent]);     

npv = 6; nc = 9; nch = 3; npf = 3; ne = length(elempart);
UDG = zeros(npv,nc,ne);
UH = zeros(nch,npf,n2(end));

filename = "run021424/soltime";

for i = 1:nproc        
    fname = filename + sprintf("%04d",time) + '_np' +num2str(i-1) +'.bin';
    fileID = fopen(fname,'r');

    ndims = fread(fileID,2,'double');                                       
    ne = nelem(i);
    UDG(:,:,elempart((n1(i)+1):(n1(i)+ne))) = reshape(fread(fileID,ndims(1),'double'),[npv nc ne]);
        
    nf = nent(i);
    tm = fread(fileID,ndims(2),'double');    
    UH(:,:,entpart((n2(i)+1):(n2(i)+nf))) = reshape(tm,[nch npf nf]);            
    
    fclose(fileID);
end
UH = reshape(UH, [nch npf*n2(end)]);

Er = UDG(:,6,:);
Ez = UDG(:,9,:);
normE = sqrt(Er.^2+Ez.^2)*3e6;

figure(1);clf;scaplot(mesh,UDG(:,1,:),[],1); axis on; colormap jet;
figure(2);clf;scaplot(mesh,UDG(:,2,:),[],1); axis on; colormap jet;
figure(3);clf;scaplot(mesh,normE,[],1); axis on; colormap jet;
figure(4);clf;scaplot(mesh,normE,[],1,1); axis on; colormap jet;

ne_star = 1e7;
param1=app.arg{1};
param2=app.arg{2};
param3=app.arg{3};
param4=app.arg{4};
charge_prefactor_tilde = param4/(param1^2*param3); 
sphi = 2e-3*ne_star*charge_prefactor_tilde*(UDG(:,2,:) - UDG(:,1,:));
alpha = 100;
sphi = sphi.*(atan(alpha*sphi)/pi + 0.5) - atan(alpha)/pi + 0.5;

l_ref = param1;
mu_ref = param2;
E_ref = param3;
De_tilde = 4.3628e-3*(normE.^0.22)/ (l_ref*mu_ref*E_ref) + sphi;

figure(5);clf;scaplot(mesh,sphi,[],1); axis on; colormap jet;
figure(6);clf;scaplot(mesh,De_tilde,[],1); axis on; colormap jet;


