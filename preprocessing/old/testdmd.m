hybrid='edg';
preconditioner = 0;
ngrid=20;
porder = 3;
nproc = 4;
elemtype=0;
nodetype=0;

% HDG
mesh = mkmesh_square(ngrid,ngrid,porder,1,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*mesh.porder);
[master,mesh] = preprocess(master,mesh,hybrid);

[epart, npart, ent2ent, ent2entStart] = meshpart(mesh,nproc,mesh.p,mesh.t);

mesh.blkSize = 2;
[mesh0,meshp0,dmd0] = domaindecomposition2(mesh,nproc);
[mesh1,meshp1,dmd1] = domaindecomposition2b(mesh,nproc);
%[mesh2,meshp2,dmd2] = domaindecomposition3(mesh,nproc);

% dmd = domaindecomposition(mesh.t2f',nproc);
%[meshp,dmd] = mkmesh_dd(mesh,nproc);
t2t = mkt2t(mesh.t,mesh.elemtype);
if strcmp(hybrid,'hdg')
    meshp = domaindecompositionmpi(mesh.t2f',t2t,mesh.t2f,nproc,preconditioner,mesh.elcon);    
else
    meshp = domaindecompositionmpi(mesh.elcon,t2t,mesh.t2f,nproc,preconditioner);
end

for i = 1:nproc
    e(1) = max(abs(meshp{i}.entpart-meshp0{i}.entpart));
    e(2) = max(abs(meshp{i}.entpartpts-meshp0{i}.entpartpts));
    e(3) = max(abs(meshp{i}.entrecv-meshp0{i}.entrecv));
    e(4) = max(abs(meshp{i}.entrecvpts-meshp0{i}.entrecvpts));
    e(5) = max(abs(meshp{i}.entsend-meshp0{i}.entsend));
    e(6) = max(abs(meshp{i}.entsendpts-meshp0{i}.entsendpts));
    e(7) = max(abs(meshp{i}.elempart-meshp0{i}.elempart));
    e(8) = max(abs(meshp{i}.elempartpts-meshp0{i}.elempartpts));
    e(9) = max(abs(meshp{i}.elemrecv-meshp0{i}.elemrecv));
    e(10) = max(abs(meshp{i}.elemrecvpts-meshp0{i}.elemrecvpts));
    e(11) = max(abs(meshp{i}.elemsend-meshp0{i}.elemsend));
    e(12) = max(abs(meshp{i}.elemsendpts-meshp0{i}.elemsendpts));
    e(13) = max(abs(meshp{i}.matrecv-meshp0{i}.matrecv));
    e(14) = max(abs(meshp{i}.matrecvpts-meshp0{i}.matrecvpts));
    e(15) = max(abs(meshp{i}.matsend-meshp0{i}.matsend));
    e(16) = max(abs(meshp{i}.matsendpts-meshp0{i}.matsendpts));
    e(17) = max(abs(meshp{i}.nbsd-meshp0{i}.nbsd));
%     e(18) = max(abs(meshp{i}.dgnodes(:)-meshp0{i}.dgnodes(:)));
%     e(19) = max(abs(meshp{i}.bf(:)-meshp0{i}.bf(:)));
    e(20) = max(abs(meshp{i}.elcon(:)-meshp0{i}.elcon(:)));
    e(21) = max(abs(meshp{i}.t2f(:)-meshp0{i}.t2f(:)));
    e(22) = max(abs(meshp{i}.cbsr_rowpts(:)-meshp0{i}.cbsr_rowpts(:)));
    e(23) = max(abs(meshp{i}.cbsr_colind(:)-meshp0{i}.cbsr_colind(:)));
    e(24) = max(abs(meshp{i}.cbsr_rowent2elem(:)-meshp0{i}.cbsr_rowent2elem(:)));
    e(25) = max(abs(meshp{i}.cbsr_colent2elem(:)-meshp0{i}.cbsr_colent2elem(:)));
    e
end

% binary file for testing
% Int ne = param[0];
% Int nf = param[1];
% Int ndof = param[2];
% Int nfe = param[3];
% Int nee = param[4];
% Int npe = param[5];
% Int npf = param[6];                
% Int preconditioner = param[11];
% Int overlappinglevel = param[12]; 
% Int hdg = param[13];
% Int nproc = param[14];

preconditioner = 0;
ngrid=9;
porder = 4;
nproc = 4;
elemtype=0;
nodetype=1;

% HDG
mesh = mkmesh_square(ngrid,ngrid,porder,1,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*mesh.porder);
[master,mesh] = preprocess(master,mesh,hybrid);
t2t = mkt2t(mesh.t,mesh.elemtype);
[epart, npart, ent2ent, ent2entStart] = meshpart(mesh,nproc,mesh.p,mesh.t);

param = zeros(20,1);
param(1) = mesh.ne;
param(2) = mesh.nf;
if strcmp(hybrid, 'hdg')
    param(3) = mesh.nf;    
    param(13) = 1;
else
    param(3) = max(mesh.elcon(:));    
    param(13) = 0;
end
param(4) = size(mesh.perm,2);
param(5) = size(t2t,2);
param(6) = size(mesh.dgnodes,1);
param(7) = size(mesh.perm,1);
param(11) = preconditioner;
param(12) = preconditioner+2;
param(14) = nproc;

fileID = fopen('dmdin.bin','w');
fwrite(fileID,param,'double');
fwrite(fileID,mesh.elcon-1,'double');
fwrite(fileID,t2t-1,'double');
fwrite(fileID,mesh.t2f-1,'double');
fwrite(fileID,epart,'double');
fwrite(fileID,npart,'double');
fclose(fileID);

% // void domaindecomposition(dmdstruct &dmd, vector<Int> &elcon, vector<Int> &elconhdg, 
% //         vector<Int> &t2t, vector<Int> &t2f, vector<Int> &elemall, vector<Int> &entall, 
% //         Int preconditioner, Int overlappinglevel, Int hdg, Int my_rank)

return;
