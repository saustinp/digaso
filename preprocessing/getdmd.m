

% if strcmp(mesh.hybrid,'hdg')                 
%     mesh.f(:,end+1) = 0;
% elseif strcmp(app.hybrid,'edg')         
%     mesh.f(:,end+1) = 1;
% elseif strcmp(app.hybrid,'iedg')    
%     isEDGface = ones(mesh.nf,1);
%     isEDGface((mesh.f(end,:)<0)) = 0;
%     mesh.f(:,end+1) = isEDGface;
% end           
%[~,~,~,facecon] = mkelcon(mesh);    
% mesh.f(:,end) = [];
facecon = 1:mesh.nf*size(mesh.perm,1);
facecon = reshape(facecon,[size(mesh.perm,1),mesh.nf]);
overlappinglevel=1;
dmd = mkdmd(mesh.t2f'-1,facecon-1,mesh.t'-1,mesh.t2f'-1,nproc,overlappinglevel,mesh.elcon-1);

hybrid = mesh.hybrid;
if strcmp(hybrid,'edg') || strcmp(hybrid,'iedg') || strcmp(hybrid,'hedg')
    mesh.blkSize = app.ncu;
elseif strcmp(hybrid,'hdg')    
    mesh.blkSize = app.ncu*master.npf;        
end
[~,meshp,dmd2] = domaindecomposition2(mesh,nproc);

err = ones(nproc,8);
for i=1:nproc    
    err(i,2)=max(abs(dmd{i}.elempart-meshp{i}.elempart));
    err(i,1)=max(abs(dmd{i}.entpart-meshp{i}.entpart));       
    err(i,3)=max(abs(dmd{i}.entsend-meshp{i}.entsend)); 
    err(i,4)=max(abs(dmd{i}.entrecv-meshp{i}.entrecv));    
    err(i,5)=max(abs(dmd{i}.elemsend-meshp{i}.elemsend));    
    err(i,6)=max(abs(dmd{i}.elemrecv-meshp{i}.elemrecv));    
    err(i,7)=max(abs(dmd{i}.matsend-meshp{i}.matsend));    
    err(i,8)=max(abs(dmd{i}.matrecv-meshp{i}.matrecv));   
    err(i,9)=max(abs(dmd{i}.elcon(:)-meshp{i}.elcon(:)));   
    err(i,10)=max(abs(dmd{i}.bcrs_rowent2elem-meshp{i}.cbsr_rowent2elem));
    err(i,11)=max(abs(dmd{i}.bcrs_colent2elem-meshp{i}.cbsr_colent2elem));
    err(i,12)=max(abs(dmd{i}.bcrs_rowent2ent-meshp{i}.cbsr_rowpts));
    err(i,13)=max(abs(dmd{i}.bcrs_colent2ent-meshp{i}.cbsr_colind));    
end

