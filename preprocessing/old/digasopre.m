function app = digasopre(app,filename,p,t,dgnodes,UDG,UH,PDG,elementtype,bndexpr,nproc,endianType,check)

if nargin < 13; check = 0; end
if nargin < 12; endianType = 0; end
if endianType == 0; endian = 'native';
elseif endianType == 1; endian = 'ieee-le';
elseif endianType == 2; endian = 'ieee-be';
end

if nargin<=2    
    % app structure
    app.elemtype = unique(elementtype);
    fileapp = [filename 'app.bin'];
    app = mkapp(app,fileapp,endian);
else   
    app.nproc = nproc;
    if isfield(app,'nfile') == 0
        app.nfile = 1;      
    end
    nfile = app.nfile;
    
    % master strcuture
    filemaster = [filename 'master.bin'];
    fileID = fopen(filemaster,'w');
    app.elemtype = unique(elementtype);
    for i = 1:length(app.elemtype)
        master = mkmasterelement(app.nd,app.porder,app.pgauss,app.pgaussR,app.elemtype(i),app.nodetype);        
        writemaster(master,fileID,endian);        
    end
    fclose(fileID);
    
    dim = app.nd;
    [t2f,t2t,f,ne,nf,nfe,nve,nvf] = mkt2ffast(t,elementtype,dim,check);
    if strcmp(app.hybrid,'hdg')         
        isEDGface = zeros(nf,1);
        hybridn = 0;
    elseif strcmp(app.hybrid,'edg')         
        isEDGface = ones(nf,1);
        hybridn = 1;
    elseif strcmp(app.hybrid,'iedg')    
        isEDGface = ones(nf,1);
        isEDGface((f(end,:)<0)) = 0;
        hybridn = 2;
    elseif strcmp(app.hybrid,'hedg')        
        isEDGface = EDGface(p,t,f);    
        hybridn = 2;    
    else
        error('app.hybrid not implemented');    
    end       

    % create elcon
    [elcon,facecon,edg] = mkelconfast(p,t,f,t2f,elementtype,isEDGface,app.porder,dim,check);        
    npfmax = size(elcon,1);    
    ndh = max(elcon(:))+1;    
    
    f = setbndnbrs(p,f,bndexpr);
    bf = reshape(f(end,t2f+1),size(t2f));
    npemax = size(dgnodes,1);
    ncd = size(dgnodes,2);
    
    if isempty(UH)==0   
        if strcmp(app.hybrid,'hdg')==1
            UH = reshape(UH,app.nch,npfmax,nf);        
        else
            UH = reshape(UH,app.nch,[]);        
        end
    end    

    ndims = zeros(30,1);
    ndims(1) = dim;  % spatial dimension
    ndims(2) = ncd; % components of dgnodes
    ndims(3) = max(nfe); % faces per element
    ndims(4) = max(nve); % vertices per element
    ndims(5) = max(nvf); % vertices per face
    ndims(6) = ne;  % number of elements
    ndims(7) = nf;  % number of faces
    ndims(8) = size(p,1);  % number of vertices
    ndims(9) = ndh;  % number of DOF of UH
    ndims(10) = npemax; % number of DG nodes per element
    ndims(11) = npfmax;  % number of DG nodes per face
    ndims(12) = app.nc;  % number of components of UDG
    ndims(13) = app.ncu; % number of compoments of U
    ndims(14) = app.ncq; % number of compoments of Q
    ndims(15) = app.ncp; % number of compoments of P
    ndims(16) = app.nch; % number of compoments of UH
    ndims(end-1) = hybridn;
    ndims(end) = nproc; % number of processors    

    mdims = zeros(10,1);
    mdims(1) = numel(UDG);
    mdims(2) = numel(UH);
    mdims(3) = numel(PDG);    
    if nproc>1                        
        elcon = reshape(elcon,[npfmax*max(nfe) ne]);
        [elem2cpu,extintelem,extintelempts,face2cpu,extintface,extintfacepts,ent2cpu,extintent,extintentpts] =...
            overlappingpartition(t2f,t,elcon,facecon,app.hybrid,app.overlappinglevel,nproc,check);
                                
        if check==1
            plotpartition(p,t+1,edg,elem2cpu,ent2cpu,extintelem,extintent,extintelempts,app.hybrid);            
            if strcmp(app.hybrid,'hdg')==1   
                dmd = mkdmd(t2f,facecon,t,t2f,nproc,app.overlappinglevel,elcon);
            else                
                dmd = mkdmd(elcon,facecon,t,t2f,nproc,app.overlappinglevel);
            end
            for i = 1:nproc
                disp(['---------------------------']);
                disp(['CPU:  ' num2str(i)]);
                disp(['---------------------------']);            
                fn = ['poioutdmd' num2str(i-1)];
                dmdcpp = readDMDSTRUCTfromBinaryFile(fn);
                validate_dmd(dmd{i},dmdcpp);            
            end
        end
                 
        if rem(nproc,nfile)~=0
            error('Number of files must be divisible by number of processors');            
        end
        nprocperfile = nproc/nfile; % number of processors per file
        meshsize = zeros(nprocperfile,nfile);
        solsize = zeros(nprocperfile,nfile);                
        app.meshsize = zeros(nprocperfile+1,nfile);
        app.solsize = zeros(nprocperfile+1,nfile);                
        for ia=1:nfile
            fileID = fopen([filename,'mesh' num2str(ia) '.bin'],'w');
            fileSO = fopen([filename,'sol' num2str(ia) '.bin'],'w');        
            for ib = 1:nprocperfile
                % ith processor
                i = (ia-1)*nprocperfile+ib;                 
                elem = extintelem{i}+1;
                face = extintface{i}+1;                
                ne = length(elem);
                nf = length(face);        
                
                if strcmp(app.hybrid,'hdg')==1                    
                    ndh = npfmax*nf;
                    nent = nf;                % number of entities    
                    nin = extintfacepts(i,1); % number of interior entities    
                    nen = extintfacepts(i,2); % number of exterior entities    
                    nbsd = unique(face2cpu(face));                    
                else
                    ent = extintent{i}+1;
                    ndh = length(ent);  
                    nent = ndh;              % number of entities
                    nin = extintentpts(i,1); % number of interior entities   
                    nen = extintentpts(i,2); % number of interior entities   
                    nbsd = unique(ent2cpu(ent));                    
                end                                
                nbsd = setdiff(nbsd,i-1);                                                
                ndims(6) = ne;  % number of elements
                ndims(7) = nf;  % number of faces   
                ndims(9) = ndh; % DOF of UH 
                ndims(17) = nent;  
                ndims(18) = extintelempts(i,1); % number of interior elements
                ndims(19) = extintfacepts(i,1); % number of interior faces                    
                ndims(20) = nin;                 % number of interior entities               
                ndims(21) = length(nbsd);        % number of neighboring cpus    
                ndims(22) = extintelempts(i,2);  % number of exterior elements
                ndims(23) = extintfacepts(i,2);  % number of exterior faces     
                ndims(24) = nen;                 % number of exterior entities     
                
                 % mesh structure              
                fwrite(fileID,ndims,'double',endian);             
                fwrite(fileID,elementtype(elem),'double',endian); 
                fwrite(fileID,t(:,elem),'double',endian);         
                fwrite(fileID,dgnodes(:,:,elem),'double',endian); 
                fwrite(fileID,elcon(:,elem),'double',endian);   
                fwrite(fileID,bf(:,elem),'double',endian);
                fwrite(fileID,t2f(:,elem),'double',endian);
                fwrite(fileID,t2t(:,elem),'double',endian);
                fwrite(fileID,f(:,face),'double',endian);                
                fwrite(fileID,isEDGface(face),'double',endian);
                fwrite(fileID,nbsd,'double',endian);
                fwrite(fileID,elem-1,'double',endian);                
                fwrite(fileID,elem2cpu(elem),'double',endian);
                fwrite(fileID,face-1,'double',endian);               
                fwrite(fileID,face2cpu(face),'double',endian);                
                N = length(ndims)+(1+nve+npemax*ncd+npfmax*nfe+nfe+nfe+nfe+2)*ne+nf*(nvf+3)+nf+nf+nf+length(nbsd);
                if strcmp(app.hybrid,'hdg')==0                    
                    fwrite(fileID,ent-1,'double',endian);   
                    fwrite(fileID,ent2cpu(ent),'double',endian);
                    N = N+2*length(ent);
                end                
                meshsize(ib,ia) = N;

                % solution        
                if isempty(UDG)==0
                    mdims(1) = app.nc*npemax*ne;
                end
                if isempty(UH)==0                                        
                    mdims(2) = app.nch*ndh;
                end
                if app.wave==1
                    if isempty(PDG)==0                
                        mdims(3) = app.ncu*npemax*ne;
                    end
                end                        
                fwrite(fileSO,mdims,'double',endian);
                M = length(mdims);
                if isempty(UDG)==0
                    fwrite(fileSO,UDG(:,:,elem),'double',endian);              
                    M = M + mdims(1);        
                end
                if isempty(UH)==0           
                    if strcmp(app.hybrid,'hdg')==0   
                        fwrite(fileSO,UH(:,ent),'double',endian);   
                    else
                        fwrite(fileSO,UH(:,:,face),'double',endian);   
                    end
                    M = M + mdims(2);
                end
                if app.wave==1
                    if isempty(PDG)==0
                        fwrite(fileSO,PDG(:,:,elem),'double',endian);
                        M = M + mdims(3);
                    end
                end                        
                solsize(ib,ia) = M;
            end
            fclose(fileID);    
            fclose(fileSO);            
            app.meshsize(:,ia) = [0; cumsum(meshsize(:,ia))];
            app.solsize(:,ia) = [0; cumsum(solsize(:,ia))];                        
        end              
        fileapp = [filename 'app.bin'];
        app = mkapp(app,fileapp,endian);                
    else    
        fileID = fopen([filename,'mesh.bin'],'w');
        fileSO = fopen([filename,'sol.bin'],'w');        
        % mesh                
        fwrite(fileID,ndims,'double',endian);            
        fwrite(fileID,elementtype,'double',endian); % local dgnodes on each subdomain        
        fwrite(fileID,t,'double',endian); % local dgnodes on each subdomain        
        fwrite(fileID,dgnodes,'double',endian); % local dgnodes on each subdomain        
        fwrite(fileID,elcon,'double',endian);
        fwrite(fileID,bf,'double',endian);
        fwrite(fileID,t2f,'double',endian);
        fwrite(fileID,t2t,'double',endian);
        fwrite(fileID,f,'double',endian);
        fwrite(fileID,isEDGface,'double',endian);
        fclose(fileID);
        meshsize = length(ndims)+(1+nve+npemax*ncd+npfmax*nfe+nfe+nfe+nfe)*ne+(nvf+3)*nf+nf;
        
        % solution
        fwrite(fileSO,mdims,'double',endian);
        solsize = length(mdims);
        if isempty(UDG)==0
            fwrite(fileSO,UDG,'double',endian);   
            solsize = solsize + numel(UDG);        
        end
        if isempty(UH)==0
            fwrite(fileSO,UH,'double',endian);     
            solsize = solsize + numel(UH);        
        end
        if app.wave==1
            if isempty(PDG)==0
                fwrite(fileSO,PDG,'double',endian);     
                solsize = solsize + numel(PDG);
            end
        end    
        fclose(fileSO);
        
        app.meshsize = meshsize;
        app.solsize = solsize;    
        fileapp = [filename 'app.bin'];
        app = mkapp(app,fileapp,endian);            
    end    
end


function app = mkapp(app,filename,endian)

if strcmp(app.appname,'euler')
    appname = 0;    
elseif strcmp(app.appname,'ns')
    appname = 1;        
elseif strcmp(app.appname,'poisson')
    appname = 2;   
elseif strcmp(app.appname,'ransSA')
    appname = 3;  
else
    error('app.appname not implemented');
end
if strcmp(app.hybrid,'hdg') 
    hybridn = 0;
elseif strcmp(app.hybrid,'edg') 
    hybridn = 1;
elseif strcmp(app.hybrid,'iedg')    
    hybridn = 2;
elseif strcmp(app.hybrid,'hedg')        
    hybridn = 2;    
else
    error('app.hybrid not implemented');    
end

if isfield(app,'wave') == 0
    app.wave = 0;
end
if isfield(app,'alag') == 0
    app.alag = 0;
end
if isfield(app,'adjoint') == 0
    app.adjoint = 0;
end
if isfield(app,'linearproblem') == 0
    app.linearproblem = 0;
end
if isfield(app,'flg_p') == 0
    app.flg_p = 0;
end
if isfield(app,'flg_g') == 0
    app.flg_g = 0;
end
if isfield(app,'print') == 0
    app.print = 0;
end
if isfield(app,'fc_u') == 0
    app.fc_u = 0;
end
if isfield(app,'fc_p') == 0
    app.fc_p = 0;
end
if isfield(app,'time') == 0
    app.time = 0;
end
if isfield(app,'dtfc') == 0
    app.dtfc = 0;
end
if isfield(app,'alpha') == 0
    app.alpha = 0;
end
if isfield(app,'timeatchingscheme') == 0
    app.timeatchingscheme = 0;
end
if isfield(app,'torder') == 0
    app.torder = 1;
end
if isfield(app,'nstage') == 0
    app.nstage = 1;
end
if isfield(app,'convStabMethod') == 0
    app.convStabMethod = 0;
end
if isfield(app,'diffStabMethod') == 0
    app.diffStabMethod = 0;
end
if isfield(app,'rotatingFrame') == 0
    app.rotatingFrame = 0;
end
if isfield(app,'viscosityModel') == 0
    app.viscosityModel = 0;
end
if isfield(app,'SGSmodel') == 0
    app.SGSmodel = 0;
end
if isfield(app,'ALE') == 0
    app.ALE = 0;
end
if isfield(app,'AV') == 0
    app.AV = 0;
end
if isfield(app,'nonlinearsolver') == 0
    app.nonlinearsolver = 0;
end
if isfield(app,'linearsolver') == 0
    app.linearsolver = 0;
end
if isfield(app,'overlappinglevel') == 0
    app.overlappinglevel = 0;
end
if isfield(app,'preconditionertype') == 0
    app.preconditionertype = 0;
end
if isfield(app,'preconditionerlevel') == 0
    app.preconditionerlevel = 1;
end
if isfield(app,'preconditionerside') == 0
    app.preconditionerside = 0;
end
if isfield(app,'newtoniter') == 0
    app.newtoniter = 10;
end
if isfield(app,'newtontol') == 0
    app.newtontol = 1e-7;
end
if isfield(app,'gmresiter') == 0
    app.gmresiter = 200;
end
if isfield(app,'restart') == 0
    app.restart = 20;
end
if isfield(app,'gmrestol') == 0
    app.gmrestol = 1e-6;
end
if isfield(app,'reuseOrdering') == 0
    app.reuseOrdering = 0;
end
if isfield(app,'reuseJacobian') == 0
    app.reuseJacobian = 0;
end
if isfield(app,'reuseResidual') == 0
    app.reuseResidual = 0;
end
if isfield(app,'jacobianStep') == 0
    app.jacobianStep = 0; 
end
if isfield(app,'orderingStep') == 0
    app.orderingStep = 0;
end
if isfield(app,'quasiNewton') == 0
    app.quasiNewton = 0;
end
if isfield(app,'quasiNewtonAccuracy') == 0
    app.quasiNewtonAccuracy = 1;  
end
if isfield(app,'orthogMethod') == 0
    app.orthogMethod = 1;
end
if isfield(app,'reorderMethod') == 0
    app.reorderMethod = 1;
end
if isfield(app,'schurImplementation') == 0
    app.schurImplementation = 1;  
end
if isfield(app,'matvecImplementation') == 0
    app.matvecImplementation = 1; 
end
if isfield(app,'precSolveImplementation') == 0
    app.precSolveImplementation = 1;  
end
if isfield(app,'precPrecision') == 0
    app.precPrecision = 1;        
end
if isfield(app,'matvecPrecision') == 0
    app.matvecPrecision = 1;      
end
if isfield(app,'orthogPrecision') == 0
    app.orthogPrecision = 1;      
end
if isfield(app,'adaptiveGMREStol') == 0
    app.adaptiveGMREStol = 0;         
end
if isfield(app,'nfile') == 0
    app.nfile = 1;      
end
if isfield(app,'flag') == 0
    app.flag = [];         
end
if isfield(app,'factor') == 0
    app.factor = [];         
end
if isfield(app,'problem') == 0
    app.problem = [];         
end
if isfield(app,'solversparam') == 0
    app.solversparam = [];         
end

app.flag   = [app.tdep app.wave app.alag app.adjoint app.linearproblem app.flg_q app.flg_p app.flg_g app.print app.flag];
app.factor = [app.fc_u app.fc_q app.fc_p app.time app.dtfc app.alpha app.factor];
app.problem = [hybridn appname app.timeatchingscheme app.torder app.nstage app.convStabMethod...
               app.diffStabMethod app.rotatingFrame app.viscosityModel app.SGSmodel app.ALE app.AV app.problem];
app.physicsparam  = reshape(cell2mat(app.arg),[],1); 
app.solversparam = [app.nonlinearsolver app.linearsolver...
                    app.overlappinglevel app.preconditionertype app.preconditionerlevel app.preconditionerside...
                    app.newtoniter app.newtontol app.gmresiter app.restart app.gmrestol...
                    app.jacobianStep app.orderingStep app.reuseOrdering app.reuseJacobian ...
                    app.reuseResidual app.quasiNewton app.quasiNewtonAccuracy app.orthogMethod ...
                    app.reorderMethod app.schurImplementation app.matvecImplementation app.precSolveImplementation...
                    app.precPrecision app.matvecPrecision app.orthogPrecision app.adaptiveGMREStol app.solversparam];
                                  
ndims = zeros(30,1);
ndims(1) = length(app.meshsize(:));
ndims(2) = length(app.solsize(:));
ndims(3) = length(app.porder(:)); 
ndims(4) = length(app.elemtype(:)); 
ndims(5) = length(app.nodetype(:)); 
ndims(6) = length(app.pgauss(:)); 
ndims(7) = length(app.pgaussR(:)); 
ndims(8) = length(app.quadtype(:)); 
ndims(9) = length(app.bcm(:)); % number of boundaries
ndims(10) = length(app.bcs(:)); % number of boundaries
ndims(11) = length(app.bcd(:)); % number of boundaries
ndims(12) = length(app.bcv(:)); % number of boundaries
ndims(13) = length(app.dt(:)); % number of time steps
ndims(14) = length(app.flag(:));  % length of flag
ndims(15) = length(app.factor(:)); % length of factor
ndims(16) = length(app.problem(:)); % length of physics
ndims(17) = length(app.physicsparam(:)); % number of physical parameters
ndims(18) = length(app.solversparam(:)); % number of solver parameters
ndims(19) = app.nproc;  % number of processors
ndims(20) = app.nfile;  % number of files to read and write (mesh, solution)
app.ndims = ndims;

fileID = fopen(filename,'w');
fwrite(fileID,app.ndims(:),'double',endian);
fwrite(fileID,app.meshsize(:),'double',endian);
fwrite(fileID,app.solsize(:),'double',endian);
fwrite(fileID,app.porder(:),'double',endian);
fwrite(fileID,app.elemtype(:),'double',endian);
fwrite(fileID,app.nodetype(:),'double',endian);
fwrite(fileID,app.pgauss(:),'double',endian);
fwrite(fileID,app.pgaussR(:),'double',endian);
fwrite(fileID,app.quadtype(:),'double',endian);
fwrite(fileID,app.bcm(:),'double',endian);
fwrite(fileID,app.bcs','double',endian);
fwrite(fileID,app.bcd(:),'double',endian);
fwrite(fileID,app.bcv','double',endian);
fwrite(fileID,app.dt(:),'double',endian);
fwrite(fileID,app.flag(:),'double',endian);
fwrite(fileID,app.factor(:),'double',endian);
fwrite(fileID,app.problem(:),'double',endian);
fwrite(fileID,app.physicsparam(:),'double',endian);
fwrite(fileID,app.solversparam(:),'double',endian);
fclose(fileID);

function writemaster(master,fileID,endian)

nfe = length(master.plocfc);
nve = length(master.permnode);
nle = size(master.permedge,2);
nplmax = size(master.permedge,1);
ndims = zeros(20,1);
ndims(1) = master.nd;
ndims(2) = master.elemtype;
ndims(3) = master.nodetype;
ndims(4) = master.npv;
ndims(5) = master.ngv;
ndims(6) = master.ngvR;
ndims(7) = nfe; 
ndims(8) = nve; 
ndims(9) = nle; 
ndims(10) = nplmax; 

% write master structure to files
fwrite(fileID,ndims(:),'double',endian);
% fwrite(fileID,master.porder(:),'double',endian);
% fwrite(fileID,master.pgauss(:),'double',endian);
% fwrite(fileID,master.pgaussR(:),'double',endian);
fwrite(fileID,master.npf(:),'double',endian);
fwrite(fileID,master.ngf(:),'double',endian);
fwrite(fileID,master.ngfR(:),'double',endian);
fwrite(fileID,master.nvf(:),'double',endian);
fwrite(fileID,master.npl(:),'double',endian);
fwrite(fileID,master.plocvl(:),'double',endian);
fwrite(fileID,master.gpvl(:),'double',endian);
fwrite(fileID,master.gwvl(:),'double',endian);
fwrite(fileID,master.shapvt(:),'double',endian);
fwrite(fileID,master.shapvg(:),'double',endian);
fwrite(fileID,master.shapvgdotshapvl(:),'double',endian);
fwrite(fileID,master.gpvlR(:),'double',endian);
fwrite(fileID,master.gwvlR(:),'double',endian);
fwrite(fileID,master.shapvtR(:),'double',endian);
fwrite(fileID,master.shapvgR(:),'double',endian);
fwrite(fileID,master.permnode(:)-1,'double',endian);
fwrite(fileID,master.permedge(:)-1,'double',endian);
%fwrite(fileID,master.philocvl(:),'double',endian);
for i = 1:nfe
    fwrite(fileID,master.plocfc{i}(:),'double',endian);
    fwrite(fileID,master.gpfc{i}(:),'double',endian);
    fwrite(fileID,master.gwfc{i}(:),'double',endian);
    fwrite(fileID,master.shapft{i}(:),'double',endian);
    fwrite(fileID,master.shapfg{i}(:),'double',endian);
    fwrite(fileID,master.shapfgdotshapfc{i}(:),'double',endian);
    fwrite(fileID,master.gpfcR{i}(:),'double',endian);
    fwrite(fileID,master.gwfcR{i}(:),'double',endian);
    fwrite(fileID,master.shapftR{i}(:),'double',endian);
    fwrite(fileID,master.shapfgR{i}(:),'double',endian);
    fwrite(fileID,master.permface{i}(:)-1,'double',endian);
    fwrite(fileID,master.face{i}(:)-1,'double',endian);
    %fwrite(fileID,master.philocfc{i}(:)-1,'double',endian);
end
