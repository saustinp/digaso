function [mesh,meshp] = writeBinaryFile(filename,mesh,master,app,UDG,UH,PDG,...
    nproc,RASlevel,DDobjectiveFlag,ent2entWeight,endianType)

if nargin<8
    nproc = 1;
end

% d: dimension
% e: element
% f: face
% v: vertex
% c: component
% p: polynomial
% g: gauss
% m: geometry
% t: time
% n: number 
% s: subdomain
% b: boundary

if length(app.dt) ~= app.ndt; error('app.ndt value does not match app.dt length.'); end

if nargin < 12; endianType = 0; end

if endianType == 0; endian = 'native';
elseif endianType == 1; endian = 'ieee-le';
elseif endianType == 2; endian = 'ieee-be';
end

if nproc == 1
    meshp = [];
    
    ndims = zeros(100,1);
    ndims(1) = mesh.nd;  % spatial dimension
    ndims(2) = mesh.ncd; % components of dgnodes
    ndims(3) = mesh.nfe; % faces per element
    ndims(4) = mesh.nve; % vertices per element
    ndims(5) = mesh.nvf; % vertices per face
    ndims(6) = mesh.ne;  % number of elements
    ndims(7) = mesh.nf;  % number of faces
    ndims(8) = mesh.nv;  % number of vertices
    ndims(9) = mesh.ndh; % degrees of UH
    ndims(10) = mesh.npe; % number of points per elements
    ndims(11) = mesh.npf; % number of points per face
    ndims(12) = mesh.nme; % geometry
    ndims(13) = mesh.nmf; % geometry
    ndims(14) = mesh.nge; % Gauss
    ndims(15) = mesh.ngf;
    ndims(16) = mesh.porder; % solution order
    ndims(17) = mesh.morder; % geometry order 
    ndims(18) = app.torder; % time order    
    ndims(19) = app.nstage; % number of DIRK stages
    ndims(20) = app.nc;
    ndims(21) = app.ncu;
    ndims(22) = app.ncq;
    ndims(23) = app.ncp;
    ndims(24) = app.nch;
    ndims(25) = app.ns; % number of subdomains
    ndims(26) = app.nb; % number of boundaries
    ndims(27) = app.ndt;
    ndims(28) = app.nparam;
    ndims(29) = app.nflag;
    ndims(30) = app.nfactor;
    ndims(31) = app.nproblem;
    ndims(32) = mesh.cbsr_nrows;
    ndims(33) = mesh.cbsr_nblks;
    ndims(34) = mesh.maxBlocksPerRow;
    ndims(35) = mesh.blkSize;    
    ndims(40) = nproc;

    fileID = fopen([filename,'.bin'],'w');

    % dimensions
    fwrite(fileID,ndims,'double',endian);
    
    % master
    fwrite(fileID,master.plocvl,'double',endian);
    fwrite(fileID,master.gpvl,'double',endian);
    fwrite(fileID,master.gwvl,'double',endian);
    fwrite(fileID,master.plocfc,'double',endian);
    fwrite(fileID,master.gpfc,'double',endian);
    fwrite(fileID,master.gwfc,'double',endian);
    fwrite(fileID,master.shapvt,'double',endian);
    fwrite(fileID,master.shapvg,'double',endian);
    fwrite(fileID,master.shapvgdotshapvl,'double',endian);
    fwrite(fileID,master.shapft,'double',endian);
    fwrite(fileID,master.shapfg,'double',endian);
    fwrite(fileID,master.shapfgdotshapfc,'double',endian);
    fwrite(fileID,master.shapmv,'double',endian);
    fwrite(fileID,master.shapmf,'double',endian);

    % app    
    fwrite(fileID,app.bcm(:),'double',endian);
    fwrite(fileID,app.bcs','double',endian);
    fwrite(fileID,app.bcd(:),'double',endian);
    fwrite(fileID,app.bcv','double',endian);
    fwrite(fileID,app.dt(:),'double',endian);
    fwrite(fileID,app.param(:),'double',endian);
    fwrite(fileID,app.flag(:),'double',endian);
    fwrite(fileID,app.factor(:),'double',endian);
    fwrite(fileID,app.problem(:),'double',endian);
    
    % mesh
    fwrite(fileID,mesh.dgnodes,'double',endian); % local dgnodes on each subdomain
    fwrite(fileID,mesh.elcon,'double',endian);
    fwrite(fileID,mesh.bf,'double',endian);
    fwrite(fileID,mesh.t2f,'double',endian);
    fwrite(fileID,mesh.perm,'double',endian);
    fwrite(fileID,mesh.permgeom,'double',endian);
    
    % solution
    fwrite(fileID,UDG,'double',endian);
    fwrite(fileID,UH,'double',endian);    
    if app.wave==1
        fwrite(fileID,PDG,'double',endian);
    end
    
    % sys
    fwrite(fileID,mesh.cbsr_rowpts,'double',endian);
    fwrite(fileID,mesh.cbsr_colind,'double',endian);

    fclose(fileID);
    
    % Compute memory requirements:
    ncu = app.ncu;
    nch = app.nch;
    ne = mesh.ne;
    nfe = mesh.nfe;
    npv = master.npe;
    npf = master.npf;
    blkSize = mesh.blkSize;
    cbsr_nblks = mesh.cbsr_nblks;

    memory = zeros(4,1);
    memory(1) = 8*blkSize*blkSize*cbsr_nblks/(1024^2);   % Hg
    memory(2) = 8*blkSize*blkSize*cbsr_nblks/(1024^2);   % Mg
    memory(3) = 8*blkSize*cbsr_nblks/(1024^2);             % v
    memory(4) = 8*npv*ncu*npf*nfe*nch*ne/(1024^2);       % DinvF
    memory(5) = 4*npv*ncu*npv*ncu*ne/(1024^2);       % Dinv
    memory(6) = 4*nch*npf*nfe*npv*ncu*ne/(1024^2);       % K
    
    disp(' ');
    disp(' ');
    disp('************************************************************');
    disp('************** Summary of memory requirements **************');
    disp('************************************************************');
    disp(' ');
    disp(['Hg: ', num2str(memory(1)), ' MB']);
    disp(['Mg: ', num2str(memory(2)), ' MB']);
    disp(['v: Restart x ', num2str(memory(3)),' MB']);
    disp(['DinvF: ', num2str(memory(4)),' MB']);
    disp(['Dinv: ', num2str(memory(5)),' MB']);
    disp(['K: ', num2str(memory(6)),' MB']);
    disp(' ');
    totalMemory1 = (memory(1) + memory(2) + memory(4)) / 1024;
    totalMemory2 = memory(3) / 1024;
    totalMemory3 = (memory(1) + memory(2) + memory(4) + memory(5) + memory(6)) / 1024;
    disp(['TOTAL MEMORY REQUIREMENTS (Newton): ', num2str(totalMemory1), ' GB  + Restart x ', num2str(totalMemory2), ' GB.']);
    disp(['TOTAL MEMORY REQUIREMENTS (quasi-Newton): ', num2str(totalMemory3), ' GB  + Restart x ', num2str(totalMemory2), ' GB.']);
    
elseif nproc > 1
    if nargin<9
        RASlevel = 0;
    end
    if nargin<10
        DDobjectiveFlag = 0;
    end
    
    if RASlevel==0
        % Either 0- or 1-face overlap
        if nargin < 11            
            [mesh,meshp] = domaindecomposition2(mesh,nproc,DDobjectiveFlag);            
        else
            [mesh,meshp] = domaindecomposition2(mesh,nproc,DDobjectiveFlag,ent2entWeight);
        end
    elseif RASlevel==1
        % 1-element overlap
        if nargin < 11
            [mesh,meshp] = domaindecomposition2b(mesh,nproc,DDobjectiveFlag);
        else
            [mesh,meshp] = domaindecomposition2b(mesh,nproc,DDobjectiveFlag,ent2entWeight);
        end
    elseif RASlevel==2
        % 2-face overlap
        if nargin < 11
            [mesh,meshp] = domaindecomposition3(mesh,nproc,DDobjectiveFlag);
        else
            [mesh,meshp] = domaindecomposition3(mesh,nproc,DDobjectiveFlag,ent2entWeight);
        end
    else
        error('RASlevel has invalid value. Only 0, 1- and 2-element overlap are currently implemented');
    end
    
%     return;
    
    nelem = zeros(1,nproc);
    memory = zeros(nproc,6);
    for i = 1:nproc
        nelem(i) = meshp{i}.ne;
        fileID = fopen([filename num2str(i-1) '.bin'],'w');
        
        ndims = zeros(100,1);
        ndims(1) = meshp{i}.nd;
        ndims(2) = meshp{i}.ncd; % components of dgnodes
        ndims(3) = meshp{i}.nfe; % faces per element
        ndims(4) = meshp{i}.nve; % vertices per element
        ndims(5) = meshp{i}.nvf; % vertices per face
        ndims(6) = meshp{i}.ne;
        ndims(7) = meshp{i}.nf;
        ndims(8) = meshp{i}.nv;
        ndims(9) = meshp{i}.ndh; % degrees of UH
        ndims(10) = master.npe; % solution
        ndims(11) = master.npf;
        ndims(12) = master.nme; % geometry
        ndims(13) = master.nmf;
        ndims(14) = master.nge; % Gauss
        ndims(15) = master.ngf;
        ndims(16) = app.porder; % solution order
        ndims(17) = app.morder; % geometry order 
        ndims(18) = app.torder; % time order    
        ndims(19) = app.nstage; % number of DIRK stages
        ndims(20) = app.nc;
        ndims(21) = app.ncu;
        ndims(22) = app.ncq;
        ndims(23) = app.ncp;
        ndims(24) = app.nch;
        ndims(25) = app.ns; % number of subdomains
        ndims(26) = app.nb; % number of boundaries
        ndims(27) = app.ndt;
        ndims(28) = app.nparam;
        ndims(29) = app.nflag;
        ndims(30) = app.nfactor;
        ndims(31) = app.nproblem;
        ndims(32) = meshp{i}.cbsr_nrows;
        ndims(33) = meshp{i}.cbsr_nblks;
        ndims(34) = meshp{i}.maxBlocksPerRow;
        ndims(35) = meshp{i}.blkSize;
        ndims(36) = meshp{i}.BJ_nrows;
        ndims(37) = meshp{i}.BJ_nblks;
        ndims(38) = meshp{i}.BK_nrows;
        ndims(39) = meshp{i}.BK_nblks;
        ndims(40) = nproc;
        ndims(41) = length(meshp{i}.entpartpts(:));
        ndims(42) = length(meshp{i}.entrecv(:));
        ndims(43) = length(meshp{i}.entsend(:));
        ndims(44) = length(meshp{i}.elempartpts(:));
        ndims(45) = length(meshp{i}.elemrecv(:));
        ndims(46) = length(meshp{i}.elemsend(:));
        ndims(47) = length(meshp{i}.nbsd(:));
        ndims(48) = length(mesh.cbsr_rowpts(:));    % This variable is no longer required in the C++ cpde
        ndims(49) = length(mesh.cbsr_colind(:));    % This variable is no longer required in the C++ code
        ndims(50) = length(meshp{i}.matrecv(:));
        ndims(51) = length(meshp{i}.matsend(:));
        
        % dimensions
        fwrite(fileID,ndims,'double',endian);
        N = 100;
        
        % master
        fwrite(fileID,master.plocvl,'double',endian);
        fwrite(fileID,master.gpvl,'double',endian);
        fwrite(fileID,master.gwvl,'double',endian);
        fwrite(fileID,master.plocfc,'double',endian);
        fwrite(fileID,master.gpfc,'double',endian);
        fwrite(fileID,master.gwfc,'double',endian);
        fwrite(fileID,master.shapvt,'double',endian);
        fwrite(fileID,master.shapvg,'double',endian);
        fwrite(fileID,master.shapvgdotshapvl,'double',endian);
        fwrite(fileID,master.shapft,'double',endian);
        fwrite(fileID,master.shapfg,'double',endian);
        fwrite(fileID,master.shapfgdotshapfc,'double',endian);
        fwrite(fileID,master.shapmv,'double',endian);
        fwrite(fileID,master.shapmf,'double',endian);
        
        N = N+length(master.plocvl(:))+length(master.gpvl(:))+length(master.gwvl(:))...
            +length(master.plocfc(:))+length(master.gpfc(:))+length(master.gwfc(:))...
            +length(master.shapvt(:))+length(master.shapvg(:))+length(master.shapvgdotshapvl(:))...
            +length(master.shapft(:))+length(master.shapfg(:))+length(master.shapfgdotshapfc(:))...
            +length(master.shapmv(:))+length(master.shapmf(:));
        
        % app        
        fwrite(fileID,app.bcm(:),'double',endian);
        fwrite(fileID,app.bcs','double',endian);
        fwrite(fileID,app.bcd(:),'double',endian);
        fwrite(fileID,app.bcv','double',endian);
        fwrite(fileID,app.dt(:),'double',endian);
        fwrite(fileID,app.param(:),'double',endian);
        fwrite(fileID,app.flag(:),'double',endian);
        fwrite(fileID,app.factor(:),'double',endian);
        fwrite(fileID,app.problem(:),'double',endian);

        N = N+length(app.bcm(:))+length(app.bcs(:))+length(app.bcd(:))+length(app.bcv(:))+length(app.dt(:))+...
            length(app.param(:))+length(app.flag(:))+length(app.factor(:))+length(app.problem(:));
        
        % mesh
        fwrite(fileID,meshp{i}.dgnodes,'double',endian); % local dgnodes on each subdomain
        fwrite(fileID,meshp{i}.elcon,'double',endian);
        fwrite(fileID,meshp{i}.bf,'double',endian);
        fwrite(fileID,meshp{i}.t2f,'double',endian);
        fwrite(fileID,meshp{i}.perm,'double',endian);
        fwrite(fileID,meshp{i}.permgeom,'double',endian);
        
        % solution
        fwrite(fileID,UDG(:,:,meshp{i}.elempart),'double',endian);   
        N = N+numel(UDG(:,:,meshp{i}.elempart));
        if strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg') || strcmp(mesh.hybrid,'hedg')
            fwrite(fileID,UH(:,meshp{i}.entpart),'double',endian);
            N = N + numel(UH(:,meshp{i}.entpart));
        elseif strcmp(mesh.hybrid,'hdg')
            UH = reshape(UH,[app.ncu master.npf mesh.nf]);
            fwrite(fileID,UH(:,:,meshp{i}.entpart),'double',endian);
            N = N + numel(UH(:,:,meshp{i}.entpart));
        end
        if app.wave==1
            fwrite(fileID,PDG(:,:,meshp{i}.elempart),'double',endian);
            N = N + numel(PDG(:,:,meshp{i}.elempart));
        end
        
        % sys
        fwrite(fileID,meshp{i}.cbsr_rowpts,'double',endian);
        fwrite(fileID,meshp{i}.cbsr_colind,'double',endian);
%         if i == 1
%              fwrite(fileID,mesh.globalEnt2entStart,'double',endian);
%              fwrite(fileID,mesh.globalEnt2ent,'double',endian);
%         end
        fwrite(fileID,meshp{i}.BJ_rowpts,'double',endian);
        fwrite(fileID,meshp{i}.BJ_colind,'double',endian);
        fwrite(fileID,meshp{i}.BK_rowpts,'double',endian);
        fwrite(fileID,meshp{i}.BK_colind,'double',endian);
        fwrite(fileID,meshp{i}.BK_rowind,'double',endian);
        fwrite(fileID,meshp{i}.entpart,'double',endian);
        fwrite(fileID,meshp{i}.entpartpts,'double',endian);
        fwrite(fileID,meshp{i}.entrecv,'double',endian);
        fwrite(fileID,meshp{i}.entrecvpts,'double',endian);
        fwrite(fileID,meshp{i}.entsend,'double',endian);
        fwrite(fileID,meshp{i}.entsendpts,'double',endian);
        fwrite(fileID,meshp{i}.elempart,'double',endian);
        fwrite(fileID,meshp{i}.elempartpts,'double',endian);
        fwrite(fileID,meshp{i}.elemrecv,'double',endian);
        fwrite(fileID,meshp{i}.elemrecvpts,'double',endian);
        fwrite(fileID,meshp{i}.elemsend,'double',endian);
        fwrite(fileID,meshp{i}.elemsendpts,'double',endian);
        fwrite(fileID,meshp{i}.nbsd,'double',endian); 
        fwrite(fileID,meshp{i}.matrecv,'double',endian);
        fwrite(fileID,meshp{i}.matrecvpts,'double',endian);
        fwrite(fileID,meshp{i}.matsend,'double',endian);
        fwrite(fileID,meshp{i}.matsendpts,'double',endian);
        
%         if i == 1
%             for j=1:nproc
%              fwrite(fileID,meshp{j}.ent2entWeightLen,'double',endian);
%             end
%         end
    
        N=N+length(meshp{i}.dgnodes(:))+length(meshp{i}.elcon(:))+length(meshp{i}.bf(:))+...
            length(meshp{i}.t2f(:))+length(meshp{i}.perm(:))+length(meshp{i}.permgeom(:))+...
            length(meshp{i}.cbsr_rowpts(:))+length(meshp{i}.cbsr_colind(:))+...
            length(meshp{i}.BJ_rowpts(:))+length(meshp{i}.BJ_colind(:))+length(meshp{i}.BK_rowpts(:))+...
            length(meshp{i}.BK_colind(:))+length(meshp{i}.BK_rowind(:))+...
            length(meshp{i}.entpart(:))+length(meshp{i}.entpartpts(:))+...
            length(meshp{i}.elempart(:))+length(meshp{i}.elempartpts(:))+...
            length(meshp{i}.entrecv(:))+length(meshp{i}.entrecvpts(:))+...
            length(meshp{i}.elemrecv(:))+length(meshp{i}.elemrecvpts(:))+...
            length(meshp{i}.entsend(:))+length(meshp{i}.entsendpts(:))+...
            length(meshp{i}.elemsend(:))+length(meshp{i}.elemsendpts(:))+length(meshp{i}.nbsd(:))+...
            length(meshp{i}.matrecv(:))+length(meshp{i}.matrecvpts(:))+...
            length(meshp{i}.matsend(:))+length(meshp{i}.matsendpts(:));
        
%         if i == 1
%             N = N + length(mesh.globalEnt2entStart) + length(mesh.globalEnt2ent) + nproc;
%         end

        % Compute memory requirements:
        ncu = app.ncu;
        nch = app.nch;
        ne = meshp{i}.ne;
        nfe = meshp{i}.nfe;
        npv = master.npe;
        npf = master.npf;
        blkSize = meshp{i}.blkSize;
        cbsr_nblks = meshp{i}.cbsr_nblks;
        
        memory(i,1) = 8*blkSize*blkSize*cbsr_nblks/(1024^2);   % Hg
        memory(i,2) = 8*blkSize*blkSize*cbsr_nblks/(1024^2);   % Mg
        memory(i,3) = 8*blkSize*cbsr_nblks/(1024^2);             % v
        memory(i,4) = 8*npv*ncu*npf*nfe*nch*ne/(1024^2);       % DinvF
        memory(i,5) = 4*npv*ncu*npv*ncu*ne/(1024^2);       % Dinv
        memory(i,6) = 4*nch*npf*nfe*npv*ncu*ne/(1024^2);       % K
        
        fclose(fileID);
    end
    
    disp(' ');
    disp(' ');
    disp('************************************************************');
    disp('************** Summary of memory requirements **************');
    disp('************************************************************');
    disp(' ');
    for i=1:nproc
        disp(['Proc. #',num2str(i),' Hg: ',num2str(memory(i,1)), ...
            ' MB   ||   Mg: ',num2str(memory(i,2)),' MB   ||   v: Restart x ', ...
            num2str(memory(i,3)),' MB   ||   DinvF: ',num2str(memory(i,4)),' MB   ||   Dinv: ',...
            num2str(memory(i,5)),' MB   ||   K: ',num2str(memory(i,6)),' MB']);
    end
    disp(' ');
    totalMemory1 = (sum(memory(:,1)) + sum(memory(:,2)) + sum(memory(:,4))) / 1024;
    totalMemory2 = sum(memory(:,3)) / 1024;
    totalMemory3 = (sum(memory(:,1)) + sum(memory(:,2)) + sum(memory(:,4)) + sum(memory(:,5)) + sum(memory(:,6))) / 1024;
    disp(['TOTAL MEMORY REQUIREMENTS (Newton): ', num2str(totalMemory1), ' GB + Restart x ', num2str(totalMemory2), ' GB.']);
    disp(['TOTAL MEMORY REQUIREMENTS (quasi-Newton): ', num2str(totalMemory3), ' GB + Restart x ', num2str(totalMemory2), ' GB.']);
else
    error('Invalid number of processors');
end

% Save structures into .mat file:
% save([filename, '.mat'],'-v7.3');

return;
