
function [mesh,meshp,dmd] = domaindecomposition2(mesh,nproc,DDobjectiveFlag,ent2entWeight)
% Domain decomposition for 0- or 1-element overlap preconditioner

if nargin<3; DDobjectiveFlag = 0; end
% % % DDobjectiveFlag:
% % %       0: Minimize edge cut
% % %       1: Minimize weighted edge-cut
% % %       2: Minimize edge-cut with weights on elements and entities for
% % %          more balanced DD

if strcmp(mesh.hybrid,'hdg')
    if nargin < 4
        [mesh,meshp] = domaindecomposition_hdg(mesh,nproc,DDobjectiveFlag);
    else
        [mesh,meshp] = domaindecomposition_hdg(mesh,nproc,DDobjectiveFlag,ent2entWeight);
    end
elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
    if nargin < 4
        [mesh,meshp] = domaindecomposition_edgv2(mesh,nproc,DDobjectiveFlag);
    else
        [mesh,meshp] = domaindecomposition_edgv2(mesh,nproc,DDobjectiveFlag,ent2entWeight);
    end
end
dmd = meshp; 
% save dmd0.mat dmd;
% save temporaryFile.mat dmd; pause

% load temporaryFile.mat;
meshp = dmd;
1;
[npf,nfe] = size(mesh.perm);
for i=1:nproc % loop over each subdomain
    disp(['Processing subdomain No. ',num2str(i)]);
    meshp{i}.perm = mesh.perm;
    meshp{i}.permgeom = mesh.permgeom;
    meshp{i}.dgnodes = mesh.dgnodes(:,:,meshp{i}.elempart);
    meshp{i}.bf = mesh.bf(:,meshp{i}.elempart);
    
    % Compute local element-to-entity connectivities
    elcon = reshape(mesh.elcon,[npf nfe mesh.ne]);
    meshp{i}.elcon = elcon(:,:,meshp{i}.elempart);
    ne = length(meshp{i}.elempart);
    if strcmp(mesh.hybrid,'hdg')
        t2f = mesh.t2f(meshp{i}.elempart,:);
        for j = 1:ne
            for k = 1:nfe
                t2f(j,k) = find(t2f(j,k) == meshp{i}.entpart);
            end
        end
        for j = 1:ne
            for k = 1:nfe
                ii = elcon(:,k,meshp{i}.elempart(j)) - (mesh.t2f(meshp{i}.elempart(j),k)-1)*npf;                
                ind = ((t2f(j,k)-1)*npf+1):(t2f(j,k)*npf);
                meshp{i}.elcon(:,k,j) = ind(ii);
            end
        end
        meshp{i}.t2f = t2f;
    elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
        entMapping = zeros(max(meshp{i}.entpart),1);
        entMapping(meshp{i}.entpart) = 1:length(meshp{i}.entpart);
        meshp{i}.elcon = entMapping(meshp{i}.elcon);
        if min(meshp{i}.elcon(:)) < 1; error('Something wrong.'); end

        facesInProcessor = mesh.t2f(meshp{i}.elempart,:);
        facesInProcessor = unique(facesInProcessor);
        facesInProcessor = facesInProcessor(facesInProcessor > 0);
        faceMapping = zeros(max(facesInProcessor),1);
        faceMapping(facesInProcessor) = 1:length(facesInProcessor);
        meshp{i}.t2f = zeros(ne, nfe);
        for j = 1:ne
            for k = 1:nfe
                if mesh.t2f(meshp{i}.elempart(j),k) > 0
                    meshp{i}.t2f(j,k) = faceMapping(mesh.t2f(meshp{i}.elempart(j),k));
                else
                    error('Something wrong.');
                end
            end
        end
        if min(meshp{i}.t2f(:)) < 1; error('Something wrong.'); end
    end
    meshp{i}.elcon = reshape(meshp{i}.elcon,[npf*nfe ne]);   
    meshp{i} = block_crs(meshp{i},mesh.hybrid);
        
    if strcmp(mesh.hybrid,'hdg')
        numTotalEntities = max(meshp{i}.t2f(:));
        numRealEntities = length(unique(meshp{i}.t2f(:)));
    elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
        numTotalEntities = max(meshp{i}.elcon(:));
        numRealEntities = length(unique(meshp{i}.elcon(:)));
    end
    
    if numRealEntities ~= numTotalEntities
        error('It looks like there are ghost (missing) entities in the mesh.');       % The construction of BJ and BK from cbsr is based on the assumption that there are no ghost entities
    end
    
%     nownent = sum(meshp{i}.entpartpts(1:2));
%     meshp{i}.BJ_rowpts = zeros(nownent+1,1);
%     meshp{i}.BJ_colind = [];
%     meshp{i}.BK_rowpts = 0;
%     meshp{i}.BK_colind = [];
%     meshp{i}.BK_rowind = 0;
%     m = 1;
%     for j = 1:nownent
%         nj = (meshp{i}.cbsr_rowpts(j)+1):meshp{i}.cbsr_rowpts(j+1);
%         cj = meshp{i}.cbsr_colind(nj);
%         i1 = find(cj<=nownent);
%         if i1(end) ~= length(i1)
%             error('ent2ent has neighboring entities in the middle of a row, instead of at the end');
%         end
%         c1 = cj(i1);
%         meshp{i}.BJ_rowpts(j+1) = meshp{i}.BJ_rowpts(j)+length(c1);
%         meshp{i}.BJ_colind = [meshp{i}.BJ_colind; c1];
%         i2 = find(cj>nownent);
%         c2 = cj(i2);
%         if ~isempty(c2)
%             meshp{i}.BK_rowpts(m+1) = meshp{i}.BK_rowpts(m)+length(c2);
%             meshp{i}.BK_colind = [meshp{i}.BK_colind; c2];
%             meshp{i}.BK_rowind(m)  = j;
%             m = m + 1;
%         end
%     end    
    nownent = sum(meshp{i}.entpartpts(1:2));
    meshp{i}.BJ_rowpts = zeros(nownent+1,1);
    meshp{i}.BJ_colind = 0*meshp{i}.cbsr_colind;
    meshp{i}.BK_rowpts = 0;
    meshp{i}.BK_colind = [];
    meshp{i}.BK_rowind = 0;
    n = 1;
    m = 1;
    for j = 1:nownent
        nj = (meshp{i}.cbsr_rowpts(j)+1):meshp{i}.cbsr_rowpts(j+1);
        cj = meshp{i}.cbsr_colind(nj);
        i1 = (cj<=nownent);
        c1 = cj(i1);
        lc1 = length(c1);
        if find(i1==0,1,'first') ~= lc1+1
            error('ent2ent has neighboring entities in the middle of a row, instead of at the end');
        end
        meshp{i}.BJ_rowpts(j+1) = meshp{i}.BJ_rowpts(j)+lc1;        
        meshp{i}.BJ_colind(n:1:(n+lc1-1)) = c1;
        n = n + lc1;
        c2 = cj(i1==0);
        if ~isempty(c2)
            meshp{i}.BK_rowpts(m+1) = meshp{i}.BK_rowpts(m)+length(c2);
            meshp{i}.BK_colind = [meshp{i}.BK_colind; c2];
            meshp{i}.BK_rowind(m)  = j;
            m = m + 1;
        end
    end
    meshp{i}.BJ_colind(meshp{i}.BJ_colind==0)=[];    
    if min(meshp{i}.BJ_colind(:)) < 1; error('Something wrong.'); end
    
    meshp{i}.nd  = mesh.nd;
    meshp{i}.ne  = size(meshp{i}.dgnodes,3);
    meshp{i}.nf  = max(meshp{i}.t2f(:));
    if meshp{i}.nf ~= length(unique(meshp{i}.t2f(:))); error('Something wrong'); end
    meshp{i}.ncd = size(meshp{i}.dgnodes,2);
    meshp{i}.nfe = size(meshp{i}.t2f,2);
    meshp{i}.nve = size(mesh.t,2);
    meshp{i}.nvf = size(mesh.tlocfc,2);
    t = mesh.t(meshp{i}.elempart,:);
    p = mesh.p(unique(t(:)),:);
    meshp{i}.nv  = size(p,1);       % TODO: This number is wrong when periodic BCs are used
    if strcmp(mesh.hybrid,'hdg')
        meshp{i}.ndh = length(meshp{i}.entpart)*npf;
    elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
        meshp{i}.ndh = length(meshp{i}.entpart);
    end
    meshp{i}.blkSize = mesh.blkSize;
    meshp{i}.BJ_nrows = length(meshp{i}.BJ_rowpts)-1;
    meshp{i}.BJ_nblks = length(meshp{i}.BJ_colind);
    meshp{i}.BK_nrows = length(meshp{i}.BK_rowpts)-1;
    meshp{i}.BK_nblks = length(meshp{i}.BK_colind);

    for j = 1:length(meshp{i}.nbsd)
        meshp{i}.entsendpts(j) = length(find(meshp{i}.entsend(:,1)==meshp{i}.nbsd(j)));
        meshp{i}.entrecvpts(j) = length(find(meshp{i}.entrecv(:,1)==meshp{i}.nbsd(j)));
        meshp{i}.elemsendpts(j) = length(find(meshp{i}.elemsend(:,1)==meshp{i}.nbsd(j)));
        meshp{i}.elemrecvpts(j) = length(find(meshp{i}.elemrecv(:,1)==meshp{i}.nbsd(j)));
    end
    meshp{i}.entsend = meshp{i}.entsend(:,2);
    meshp{i}.entrecv = meshp{i}.entrecv(:,2);
    meshp{i}.elemsend = meshp{i}.elemsend(:,2);
    meshp{i}.elemrecv = meshp{i}.elemrecv(:,2);
    meshp{i}.ent2entWeightLen = 2*(meshp{i}.BJ_nrows + meshp{i}.BJ_nblks + meshp{i}.BK_nblks);
end

if strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
%     for i=1:nproc % loop over each subdomain
%         meshp{i}.matrecv = [];
%         meshp{i}.matsend = [];
%     end    
%     for i=1:nproc
%         disp(['Processing subdomain No. ',num2str(i),' (step 3/3)']);              
%         for j = 1:size(dmd{i}.entrecv,1)
%             ncpu = dmd{i}.entrecv(j,1);
%             edgj = dmd{i}.entrecv(j,2);
%             edgjg = meshp{i}.entpart(edgj);            
%             [~,em] = find(mesh.elcon(:,dmd{ncpu}.elempart)==edgjg);            
%             entj{j} = unique(dmd{ncpu}.elempart(em));
%         end
%         for j = 1:size(dmd{i}.entrecv,1)
%             ncpu = dmd{i}.entrecv(j,1);
%             edgj = dmd{i}.entrecv(j,2);
%             jr = (meshp{i}.cbsr_rowpts(edgj)+1):meshp{i}.cbsr_rowpts(edgj+1);
%             edgk = meshp{i}.cbsr_colind(jr);    
%             edgjg = meshp{i}.entpart(edgj);
%             edgkg = meshp{i}.entpart(edgk);
%             m = find(dmd{ncpu}.entsend(:,3) == edgjg); 
%             if dmd{ncpu}.entsend(m,1)~=i
%                 error('something wrong');                        
%             end
%             nedgj = dmd{ncpu}.entsend(m,2);     
%             njr = (meshp{ncpu}.cbsr_rowpts(nedgj)+1):meshp{ncpu}.cbsr_rowpts(nedgj+1);
%             nedgk = meshp{ncpu}.cbsr_colind(njr);    
%             nedgkg = meshp{ncpu}.entpart(nedgk);
%             %[~,je] = find(mesh.elcon==edgjg);                                  
%             je = entj{j};
%             for k = 1:length(edgkg)
%                 %[~,ke] = find(mesh.elcon==edgkg(k));                
%                 kj = find(dmd{i}.entrecv(:,2)==edgk(k));
%                 if isempty(kj)==0                    
%                     ke = entj{kj};
%                     %jke = intersect(je,ke);
%                     %jke = je(ismember(je,ke));    
%                     jke = je(ismemberBuiltinTypes(je,ke));                    
%                     if all(ismemberBuiltinTypes(jke,meshp{i}.elempart))==0                        
%                         meshp{i}.matrecv = [meshp{i}.matrecv; ncpu jr(k)-1];                                                                                                                            
%                         nk = find(nedgkg==edgkg(k));
%                         meshp{ncpu}.matsend = [meshp{ncpu}.matsend; i njr(nk)-1]; 
%                     end
%                 end
%             end            
%         end                    
%     end
    
    for i=1:nproc % loop over each subdomain
        nentrecv = size(dmd{i}.entrecv,1);
        nentsend = size(dmd{i}.entsend,1);
        m = 0;
        for j = 1:nentrecv
            edgj = dmd{i}.entrecv(j,2);
            jr = (meshp{i}.cbsr_rowpts(edgj)+1):meshp{i}.cbsr_rowpts(edgj+1);
            m = max([m length(jr)]);
        end
        meshp{i}.matrecv = zeros(nentrecv*m,2);
        meshp{i}.matsend = zeros(nentsend*m,2);        
    end    
    ksend = ones(nproc,1);
    krecv = ones(nproc,1);
    for i=1:nproc
        disp(['Processing subdomain No. ',num2str(i),' (step 3/3)']);          
        for j = 1:size(dmd{i}.entrecv,1)
            %ncpu = dmd{i}.entrecv(j,1);
            %edgj = dmd{i}.entrecv(j,2);
            %edgjg = meshp{i}.entpart(edgj);                             
            %[~,em] = find(mesh.elcon(:,dmd{ncpu}.elempart)==edgjg);            
            %entj{j} = unique(dmd{ncpu}.elempart(em));            
            nj = dmd{i}.entrecv(j,3);
            rj = (mesh.cbsr_rowent2elem(nj)+1):mesh.cbsr_rowent2elem(nj+1);
            entj{j} = mesh.cbsr_colent2elem(rj); % list of neirghboring elements                
        end
        for j = 1:size(dmd{i}.entrecv,1)
            ncpu = dmd{i}.entrecv(j,1);
            edgj = dmd{i}.entrecv(j,2);
            jr = (meshp{i}.cbsr_rowpts(edgj)+1):meshp{i}.cbsr_rowpts(edgj+1);
            edgk = meshp{i}.cbsr_colind(jr);    
            edgjg = meshp{i}.entpart(edgj);
            edgkg = meshp{i}.entpart(edgk);
            m = find(dmd{ncpu}.entsend(:,3) == edgjg); 
            if dmd{ncpu}.entsend(m,1)~=i; error('Something wrong.'); end
            nedgj = dmd{ncpu}.entsend(m,2);     
            njr = (meshp{ncpu}.cbsr_rowpts(nedgj)+1):meshp{ncpu}.cbsr_rowpts(nedgj+1);
            nedgk = meshp{ncpu}.cbsr_colind(njr);    
            nedgkg = meshp{ncpu}.entpart(nedgk);                                        
            je = entj{j}; % elements contain edgjg
            
            for k = 1:length(edgkg)
                kj = find(dmd{i}.entrecv(:,2)==edgk(k));
                
                % Communication of the block is only required if the other
                % column entity is not in the processor:
                if isempty(kj)==0                    
                    ke = entj{kj}; % elements contain edgkg
                    jke = je(ismemberBuiltinTypes(je,ke));  % elements contain both edgjg and edgkg
                    
                    % Communication of the block is only required if there
                    % are elements simultaneously neighboring both entities
                    % that are not in the processor:
                    % TODO: Why not if all(ismemberBuiltinTypes(jke,meshp{i}.elempart(:)))==0 ?
                    if meshp{i}.elempartpts(1) > 0
                        if all(ismemberBuiltinTypes(jke,meshp{i}.elempart(meshp{i}.elempartpts(1):end)))==0
                            meshp{i}.matrecv(krecv(i),:) = [ncpu jr(k)-1];
                            krecv(i) = krecv(i)+1;
                            nk = nedgkg==edgkg(k);
                            if length(find(nk==1))~=1; error('Something wrong.'); end
                            meshp{ncpu}.matsend(ksend(ncpu),:) = [i njr(nk)-1];
                            ksend(ncpu) = ksend(ncpu)+1;
                        end
                    end
                end
            end
        end
    end
    for i=1:nproc
        ind = meshp{i}.matrecv(:,1)==0;
        meshp{i}.matrecv(ind,:)=[];
        ind = meshp{i}.matsend(:,1)==0;
        meshp{i}.matsend(ind,:)=[];    
    end
elseif strcmp(mesh.hybrid,'hdg')
    % TODO: May not provide the right preconditioner for some meshes. Need
    % to use an approach like the one for edg and iedg to make sure that
    % does not happen.
    for i=1:nproc
        meshp{i}.matsend(:,2) = meshp{i}.cbsr_rowpts(meshp{i}.matsend(:,2));
        meshp{i}.matrecv(:,2) = meshp{i}.cbsr_rowpts(meshp{i}.matrecv(:,2));
    end    
end

for i=1:nproc
    for j = 1:length(meshp{i}.nbsd)
        meshp{i}.matsendpts(j) = length(find(meshp{i}.matsend(:,1)==meshp{i}.nbsd(j)));
        meshp{i}.matrecvpts(j) = length(find(meshp{i}.matrecv(:,1)==meshp{i}.nbsd(j)));
    end
    meshp{i}.matsend = meshp{i}.matsend(:,2);
    meshp{i}.matrecv = meshp{i}.matrecv(:,2);    
end


function [mesh,dmd] = domaindecomposition_hdg(mesh,nproc,DDobjectiveFlag,ent2entWeight)

if nargin < 3; DDobjectiveFlag = 0; end

if DDobjectiveFlag == 0
    [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart(mesh,nproc);
elseif DDobjectiveFlag == 1
    if nargin < 4;
        [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart_v2(mesh,nproc);
    else
        [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart_v2(mesh,nproc,ent2entWeight);
    end
elseif DDobjectiveFlag == 2
    [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart, ...
        mesh.elemWeightInProc, mesh.entWeightInProc] = meshpart_v3(mesh,nproc);
else
    error('DDobjectiveFlag has invalid value');
end

% load('../problem/ns/LES/T106C_HOW4/p2_BaselineMesh_CenaeroBCs/T106C_P2D2Medium_dmd.mat');
% 1;
% % plot domain partition
%plotpart(mesh, nproc, elempart, entpart);

% save('temporaryFile.mat','elempart','entpart','mesh','-v7.3'); pause
% load('temporaryFile.mat');

elempart = elempart+1;
entpart = entpart+1;

% loop over each subdomain
checkElements = 0;
checkEntities = 0;
for i=1:nproc    
    % list of elements in subdomain i                       
    indelem{i} = find(elempart==i); 
    checkElements = checkElements + length(indelem{i});
    
    % list of edg nodes in subdomain i                       
    indent{i} = find(entpart==i); 
    checkEntities = checkEntities + length(indent{i});
end

numElem = mesh.ne;
numEntities = max(mesh.t2f(:));
if numEntities ~= length(unique(mesh.t2f(:))); error('It looks like there are ghost entities in the mesh.'); end
if checkElements ~= numElem; error('Domain decomposition for elements was not performed properly.'); end
if checkEntities ~= numEntities; error('Domain decomposition for entities was not performed properly.'); end

% Check that all elements in the processor have at least one entity in the
% processor
disp(' ');
disp('Check that all elements in the processor have at least one entity in the processor...');
disp(' ');
for i=1:nproc
    disp(['Proc. No. ', num2str(i),'/',num2str(nproc)]);
    elementsToRemove = [];
    for j=1:length(indelem{i})
        if length(setdiff(indent{i},mesh.t2f(indelem{i}(j),:))) == length(indent{i})
            entitiesInElement = mesh.t2f(indelem{i}(j),:);
            candidateProcessors = entpart(entitiesInElement);
            k = mode(candidateProcessors);
            indelem{k} = sort([indelem{k};indelem{i}(j)]);
            elementsToRemove = [elementsToRemove;indelem{i}(j)];
            warning('An element in the processor had did not have any entity in the processors. The issue has been fixed successfully.');
        end
    end
    if ~isempty(elementsToRemove)
        indelem{i} = sort(setdiff(indelem{i},elementsToRemove));
    end
end

% t2t = mkt2t(mesh.t,mesh.elemtype);

% loop over each subdomain
for i=1:nproc    
    disp(['Proc. No. ', num2str(i),'/',num2str(nproc)]);
    % Find all elements connected to entities in processor (i.e. elements
    % required to compute the rows of entities in the processor)
    elemNeighboringEntInProc = [];
    for j = 1:length(indent{i})
        nj = indent{i}(j);
        rj = (mesh.cbsr_rowent2elem(nj)+1):mesh.cbsr_rowent2elem(nj+1);
        newElem = mesh.cbsr_colent2elem(rj);
        elemNeighboringEntInProc = [elemNeighboringEntInProc; newElem(:)]; 
    end
    elemNeighboringEntInProc = unique(elemNeighboringEntInProc(:));
    
    % find all interface faces
    f2f = mesh.f2f(:,indent{i});
    i1 = find(f2f(:)>0);
    in = ones(size(f2f));
    in(i1) = ismember(f2f(i1),indent{i});
    [~,j1] = find(in==0);
    j1 = unique(j1);
    
    % list of all interface faces
    intf = indent{i}(j1);
    
    % Find all elements connected to interface entities (i.e. elements
    % required to compute the rows of interface entities)
    elemNeighboringInterfEnt = [];
    for j = 1:length(intf)
        nj = intf(j);
        rj = (mesh.cbsr_rowent2elem(nj)+1):mesh.cbsr_rowent2elem(nj+1);
        newElem = mesh.cbsr_colent2elem(rj);
        elemNeighboringInterfEnt = [elemNeighboringInterfEnt; newElem(:)]; 
    end
    elemNeighboringInterfEnt = unique(elemNeighboringInterfEnt(:));
    
    % find all faces that are connected to interface faces
    f2f = f2f(:,j1);
    i1 = find(f2f(:)>0);
    in = ones(size(f2f));
    in(i1) = ismember(f2f(i1),indent{i});
    f0 = unique(f2f(in==0));  % faces do not belong to subdomain i
    %f1 = unique(f2f(in==1));
    %f1 = f1(f1>0);            % faces belong to subdomain i        
    
    % [interior faces, own interface faces, other faces]   
    dmd{i}.entpart = [setdiff(indent{i},intf); intf; f0];
    dmd{i}.entpartpts = [length(indent{i})-length(intf) length(intf) length(f0)];
        
    % store faces received from neighboring subdomains to perform the matrix-vector product
    dmd{i}.entrecv = [0*f0 length(indent{i})+(1:1:length(f0))' f0];
    for k = 1:nproc
        if k ~= i
            in = ismember(f0,indent{k});
            dmd{i}.entrecv(in,1) = k;
        end
    end
    dmd{i}.entrecv = unique(dmd{i}.entrecv,'rows');    
    dmd{i}.nbsd = unique(dmd{i}.entrecv(:,1))';        
    if dmd{i}.nbsd(1)==0
        error('something wrong');
    end
    
    % find all elements that are connected to interface faces
    elem = mesh.f(intf,end-1:end);
    i1 = find(elem(:)>0);
    in = ones(size(elem));
    in(i1) = ismember(elem(i1),indelem{i});
    elem0 = unique(elem(in==0)); % elements do not belong to subdomain i
    elem1 = unique(elem(in==1));
    elem1 = elem1(elem1>0);      % elements belong to subdomain i
    
    % [interior elements, interface elements, other elements]  
%     elemi  = setdiff(indelem{i},elem1);
%     elemj = t2t(elem1,:);    
%     elemj = unique(elemj(:));
%     i1 = elemj==0;
%     elemj(i1) = [];
%     in = ismember(elemj,elemi);
%     elem2 = unique(elemj(in==1));     
%     dmd{i}.elempart = [setdiff(elemi,elem2); [elem2; elem1]; elem0];
%     dmd{i}.elempartpts = [length(elemi)-length(elem2) length(elem2)+length(elem1) length(elem0)];    
    elem2 = setdiff(elemNeighboringInterfEnt, [elem1; elem0]);
    if ~all(ismember(elem2,indelem{i})); error('Something wrong.'); end
    dmd{i}.elempart = [setdiff(indelem{i},[elem2; elem1]); [elem2; elem1]; elem0];
    dmd{i}.elempartpts = [length(setdiff(indelem{i},[elem2; elem1])) length(elem2)+length(elem1) length(elem0)]; 
    
    % store elements received from neighboring subdomains to perform the matrix-vector assembly
    dmd{i}.elemrecv = [0*elem0 length(indelem{i})+(1:1:length(elem0))' elem0];
    for k = 1:nproc
        if k ~= i
            in = ismember(elem0,indelem{k});
            dmd{i}.elemrecv(in,1) = k;
        end
    end
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');
    
    % Check that no new boundary subdomain appear for the elements compared
    % to the entities:
    newNeighSubdomains = setdiff(dmd{i}.elemrecv(:,1), dmd{i}.nbsd);
    if ~isempty(newNeighSubdomains); error('New neighboring subdomains appear for the elements compated to the entities. The DD will not work in this case'); end
    % In case this happens, it sufficies to append the new neighboring
    % subdomains to dmd{i}.nbsd (and, maybe, sort the new vector)
    
    % Check that all elements neighboring entities in processor have been
    % detected:
    elemNotDetected = setdiff(elemNeighboringEntInProc, dmd{i}.elempart);
    if ~isempty(elemNotDetected); error('Element neighboring entity in processor have not been detected.'); end
end

% store faces sent to neighboring subdomains to perform the matrix-vector product
for k = 1:nproc
    dmd{k}.entsend = [];
end
for i = 1:nproc          
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.entrecv(:,1)==k;
        tm = dmd{i}.entrecv(ii,:);
        tm(:,1) = i;                
        for m = 1:size(tm,1)            
            tm(m,2) = find(dmd{k}.entpart==tm(m,3));
        end
        dmd{k}.entsend = [dmd{k}.entsend; tm];        
    end    
end

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc
    dmd{k}.elemsend = [];
end
for i = 1:nproc        
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.elemrecv(:,1)==k;
        tm = dmd{i}.elemrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)
            tm(m,2) = find(dmd{k}.elempart==tm(m,3));
        end
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end    
end

% store faces sent to neigboring cpus to assemble the preconditioner
for k = 1:nproc
    dmd{k}.matrecv = dmd{k}.entrecv;
    on = [];
    for i = 1:size(dmd{k}.matrecv,1)
        fi = dmd{k}.matrecv(i,3);
        ei = mesh.f(fi,end-1:end);
        if (ei(2)<=0) || (all(ismember(ei,dmd{k}.elempart))==1)
            on = [on i];                    
        end
    end    
    dmd{k}.matrecv(on,:) = [];
end
for k = 1:nproc
    dmd{k}.matsend = [];
end
for i = 1:nproc          
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.matrecv(:,1)==k;
        tm = dmd{i}.matrecv(ii,:);
        tm(:,1) = i;                
        for m = 1:size(tm,1)            
            tm(m,2) = find(dmd{k}.entpart==tm(m,3));
        end
        dmd{k}.matsend = [dmd{k}.matsend; tm];        
    end    
end

function [mesh,dmd] = domaindecomposition_edg(mesh,nproc,DDobjectiveFlag,ent2entWeight)

if nargin < 3; DDobjectiveFlag = 0; end

if DDobjectiveFlag == 0
    [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart(mesh,nproc);
elseif DDobjectiveFlag == 1
    if nargin < 4;
        [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart_v2(mesh,nproc);
    else
        [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart_v2(mesh,nproc,ent2entWeight);
    end
elseif DDobjectiveFlag == 2
    [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart, ...
        mesh.elemWeightInProc, mesh.entWeightInProc] = meshpart_v3(mesh,nproc);
else
    error('DDobjectiveFlag has invalid value');
end

% plot domain partition
% plotpart(mesh, nproc, elempart, entpart);

elempart = elempart+1;
entpart = entpart+1;

% loop over each subdomain
for i=1:nproc
    
    % list of elements in subdomain i                       
    indelem{i} = find(elempart==i); 
    
    % list of edg nodes in subdomain i                       
    indent{i} = find(entpart==i); 
end

% Check that all elements in the processor have at least one entity in the
% processor
for i=1:nproc
    elementsToRemove = [];
    for j=1:length(indelem{i})
        if length(setdiff(indent{i},mesh.elcon(:,indelem{i}(j)))) == length(indent{i})
            entitiesInElement = mesh.elcon(:,indelem{i}(j));
            candidateProcessors = entpart(entitiesInElement);
            k = mode(candidateProcessors);
            indelem{k} = sort([indelem{k};indelem{i}(j)]);
            elementsToRemove = [elementsToRemove;indelem{i}(j)];
            warning('An element in the processor had did not have any entity in the processors. The issue has been fixed successfully.');
        end
    end
    if ~isempty(elementsToRemove)
        indelem{i} = sort(setdiff(indelem{i},elementsToRemove));
    end
end

t2t = mkt2t(mesh.t,mesh.elemtype);

% loop over each subdomain
for i=1:nproc
    % find all interface edgnodes
    edgnumber = zeros(length(indent{i}),1); 
    for j = 1:length(indent{i})
        [~,j1] = find(mesh.elcon==indent{i}(j));
        j1 = unique(j1);
        edg1 = mesh.elcon(:,j1);
        edg1 = unique(edg1(:));     % list of neirghboring edg nodes
        if all(ismember(edg1,indent{i}))
            % edgnode indent{i}(j) is fully inside the subdomain i (all
            % neighboring edg nodes are in the subdomain i)
            edgnumber(j) = 2;
        else
            % edgnode indent{i}(j) is on the interface between two subdomains
            edgnumber(j) = 1;            
        end
    end
        
    % list of all interface edgnodes in the subdomain
    intf = indent{i}(edgnumber==1);
    
    % find all neighboring edgnodes that are connected to interface edgnodes
    edg0 = [];
    for j = 1:length(intf)
        [~,j1]=find(mesh.elcon==intf(j));
        j1 = unique(j1);
        edg1 = mesh.elcon(:,j1);
        edg1 = unique(edg1(:));
        in = ismember(edg1,indent{i});
        edg0 = [edg0; edg1(in==0)];
    end
    edg0 = unique(edg0);        % list of all edg nodes that are not in the subdomain but are connected to edg nodes in the subdomain
    
    % [interior edgnodes, own interface edgnodes, other edgnodes (neighbors in other processors)]
    dmd{i}.entpart = [setdiff(indent{i},intf); intf; edg0];
    dmd{i}.entpartpts = [length(setdiff(indent{i},intf)) length(intf) length(edg0)];    
    
    % store edgnodes received from neighboring subdomains to perform the matrix-vector product
    dmd{i}.entrecv = [0*edg0 length(indent{i})+(1:1:length(edg0))' edg0];
    for k = 1:nproc
        if k ~= i
            in = ismember(edg0,indent{k});
            dmd{i}.entrecv(in,1) = k;
        end
    end
    dmd{i}.entrecv = unique(dmd{i}.entrecv,'rows');    
    dmd{i}.nbsd = unique(dmd{i}.entrecv(:,1))';     % list of subdomains from which info is required for matrix-vector product
    
    % find all elements that are connected to interface edgnodes
    elem0 = [];     % elements not in the process that neighbor edg nodes in the process
    elem1 = [];     % elements in the process that neighbor edg nodes NOT in the process
    for j = 1:length(intf)
        [~,j1]=find(mesh.elcon==intf(j));
        j1 = unique(j1);        
        in = ismember(j1,indelem{i});
        elem0 = [elem0; j1(in==0)];
        elem1 = [elem1; j1(in==1)];
    end
    elem0 = unique(elem0);
    elem1 = unique(elem1);
    in = zeros(length(elem1),1);        
    for j = 1:length(elem1)        
        if all(ismember(mesh.elcon(:,elem1(j)),indent{i}))            
            in(j) = 1;
        end
    end  
    elem1 = elem1(in==0);    
    
    % [interior elements, interface elements, other elements (neighbors in other processors)]  
    elemi  = setdiff(indelem{i},elem1);
    elemj = t2t(elem1,:);    
    elemj = unique(elemj(:));
    i1 = elemj==0;
    elemj(i1) = [];
    in = ismember(elemj,elemi);
    elem2 = unique(elemj(in==1));         
    dmd{i}.elempart = [setdiff(elemi,elem2); [elem2; elem1]; elem0];
    dmd{i}.elempartpts = [length(elemi)-length(elem2) length(elem2)+length(elem1) length(elem0)];        
%     dmd{i}.elempart = [setdiff(indelem{i},elem1); elem1; elem0];
%     dmd{i}.elempartpts = [length(setdiff(indelem{i},elem1)) length(elem1) length(elem0)];    
    
    % store elements received from neighboring subdomains to assemble linear system
    dmd{i}.elemrecv = [0*elem0 length(indelem{i})+(1:1:length(elem0))' elem0];
    for k = 1:nproc
        if k ~= i
            in = ismember(elem0,indelem{k});
            dmd{i}.elemrecv(in,1) = k;
        end
    end
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');        
end

% store edgnodes sent to neighboring subdomains to perform the matrix-vector product
for k = 1:nproc
    dmd{k}.entsend = [];
end
for i = 1:nproc
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.entrecv(:,1)==k;
        tm = dmd{i}.entrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)
            tm(m,2) = find(dmd{k}.entpart==tm(m,3));
        end
        dmd{k}.entsend = [dmd{k}.entsend; tm];        
    end    
end

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc
    dmd{k}.elemsend = [];
end
for i = 1:nproc        
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.elemrecv(:,1)==k;
        tm = dmd{i}.elemrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)            
            tm(m,2) = find(dmd{k}.elempart==tm(m,3));
        end
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end    
end

% for i = 1:nproc
%     dmd{i}.elempart'
%     dmd{i}.elemsend
%     dmd{i}.elemrecv 
%     dmd{i}.entpart'
%     dmd{i}.entsend
%     dmd{i}.entrecv
% end

function [mesh,dmd] = domaindecomposition_edgv2(mesh,nproc,DDobjectiveFlag,ent2entWeight)

if nargin < 3; DDobjectiveFlag = 0; end

if DDobjectiveFlag == 0
    [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart(mesh,nproc);
elseif DDobjectiveFlag == 1
    if nargin < 4;
        [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart_v2(mesh,nproc);
    else
        [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart_v2(mesh,nproc,ent2entWeight);
    end
elseif DDobjectiveFlag == 2
    [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart, ...
        mesh.elemWeightInProc, mesh.entWeightInProc] = meshpart_v3(mesh,nproc);
else
    error('DDobjectiveFlag has invalid value');
end

% plot domain partition
% plotpart(mesh, nproc, elempart, entpart);

% globalEnt2ent = mesh.globalEnt2ent;
% globalEnt2entStart = mesh.globalEnt2entStart;
% save('tmp_dmd_512np.mat','elempart','entpart','globalEnt2ent','globalEnt2entStart','-v7.3'); pause
% load('tmp_dmd_512np.mat');
% mesh.globalEnt2ent = globalEnt2ent;
% mesh.globalEnt2entStart = globalEnt2entStart;

elempart = elempart+1;
entpart = entpart+1;

% loop over each subdomain
checkElements = 0;
checkEntities = 0;
for i=1:nproc
    % list of elements in subdomain i                       
    indelem{i} = find(elempart==i); 
    checkElements = checkElements + length(indelem{i});
    
    % list of edg nodes in subdomain i                       
    indent{i} = find(entpart==i); 
    checkEntities = checkEntities + length(indent{i});
end

numElem = mesh.ne;
numEntities = max(mesh.elcon(:));
if numEntities ~= length(unique(mesh.elcon(:))); error('It looks like there are ghost entities in the mesh.'); end
if checkElements ~= numElem; error('Domain decomposition for elements was not performed properly.'); end
if checkEntities ~= numEntities; error('Domain decomposition for entities was not performed properly.'); end

% Check that all elements in the processor have at least one entity in the
% processor
disp(' ');
disp('Check that all elements in the processor have at least one entity in the processor...');
disp(' ');
for i=1:nproc
    disp(['Proc No. ', num2str(i), ' / ', num2str(nproc)]);
    elementsToRemove = [];
    for j=1:length(indelem{i})
        if length(setdiff(indent{i},mesh.elcon(:,indelem{i}(j)))) == length(indent{i})
            entitiesInElement = mesh.elcon(:,indelem{i}(j));
            candidateProcessors = entpart(entitiesInElement);
            k = mode(candidateProcessors);
            indelem{k} = sort([indelem{k};indelem{i}(j)]);
            elementsToRemove = [elementsToRemove;indelem{i}(j)];
            warning('An element in the processor had did not have any entity in the processors. The issue has been fixed successfully.');
        end
    end
    if ~isempty(elementsToRemove)
        indelem{i} = sort(setdiff(indelem{i},elementsToRemove));
    end
end

% t2t = mkt2t(mesh.t,mesh.elemtype);

% loop over each subdomain
for i=1:nproc    
    disp(['Preprocessing subdomain No. ',num2str(i)]);
    
    % Find all elements connected to entities in processor (i.e. elements
    % required to compute the rows of entities in the processor)
    elemNeighboringEntInProc = [];
    for j = 1:length(indent{i})
        nj = indent{i}(j);
        rj = (mesh.cbsr_rowent2elem(nj)+1):mesh.cbsr_rowent2elem(nj+1);
        newElem = mesh.cbsr_colent2elem(rj);
        elemNeighboringEntInProc = [elemNeighboringEntInProc; newElem(:)]; 
    end
    elemNeighboringEntInProc = unique(elemNeighboringEntInProc(:));
    
    % find all interface nodes
    nent = length(indent{i});
    edgnumber = zeros(length(indent{i}),1);
    for j = 1:nent
        nj = indent{i}(j);
        rj = (mesh.cbsr_rowpts(nj)+1):mesh.cbsr_rowpts(nj+1);
        edg1 = mesh.cbsr_colind(rj); % list of neirghboring edg nodes
        if all(ismember(edg1,indent{i}))
            % edgnode indent{i}(j) is fully inside the subdomain i (all
            % neighboring edg nodes are in the subdomain i)
            edgnumber(j) = 2;
        else
            % edgnode indent{i}(j) is on the interface between two subdomains
            edgnumber(j) = 1;            
        end        
    end
    
    % list of all interface edgnodes in the subdomain
    intf = indent{i}(edgnumber==1);

    % Find all elements connected to interface entities (i.e. elements
    % required to compute the rows of interface)
    elemNeighboringInterfEnt = [];
    for j = 1:length(intf)
        nj = intf(j);
        rj = (mesh.cbsr_rowent2elem(nj)+1):mesh.cbsr_rowent2elem(nj+1);
        newElem = mesh.cbsr_colent2elem(rj);
        elemNeighboringInterfEnt = [elemNeighboringInterfEnt; newElem(:)]; 
    end
    elemNeighboringInterfEnt = unique(elemNeighboringInterfEnt(:));
    
    % find all neighboring edgnodes that are connected to interface edgnodes
    edg0 = [];
    for j = 1:length(intf)
        nj = intf(j);
        rj = (mesh.cbsr_rowpts(nj)+1):mesh.cbsr_rowpts(nj+1);
        edg1 = mesh.cbsr_colind(rj); % list of neirghboring edg nodes
        edg1 = unique(edg1(:));
        in = ismember(edg1,indent{i});
        edg0 = [edg0; edg1(in==0)];                
    end
    edg0 = unique(edg0); % list of all edg nodes that are not in the subdomain but are connected to edg nodes in the subdomain
    
    % [interior edgnodes, own interface edgnodes, other edgnodes (neighbors in other processors)]
    dmd{i}.entpart = [setdiff(indent{i},intf); intf; edg0];
    dmd{i}.entpartpts = [length(setdiff(indent{i},intf)) length(intf) length(edg0)];    
    
    % store edgnodes received from neighboring subdomains to perform the matrix-vector product
    dmd{i}.entrecv = [0*edg0 length(indent{i})+(1:1:length(edg0))' edg0];
    for k = 1:nproc
        if k ~= i
            in = ismember(edg0,indent{k});
            dmd{i}.entrecv(in,1) = k;
        end
    end
    if any(dmd{i}.entrecv(:,1) == 0); error('Domain from which entity needs to be communcated was not found.'); end
    dmd{i}.entrecv = unique(dmd{i}.entrecv,'rows');       % This is required to sort by entries in first column  
    dmd{i}.nbsd = unique(dmd{i}.entrecv(:,1))';     % list of subdomains from which info is required for matrix-vector product    
    
    
%     % find a subset of elements containing interface edgnodes    
%     e = indelem{i};
%     for j = 1:length(dmd{i}.nbsd)
%         k = dmd{i}.nbsd(j);
%         e = [e; indelem{k}];
%     end
%     e = unique(e(:));
%     elc = mesh.elcon(:,e);
%     
%     % find all elements that are connected to interface edgnodes
%     elem0 = [];     % elements not in the process that neighbor edg nodes in the process
%     elem1 = [];     % elements in the process that neighbor edg nodes NOT in the process
%     for j = 1:length(intf)        
% %         [~,j1]=find(mesh.elcon==intf(j));
% %         j1 = unique(j1);        
%         [~,j2]=find(elc==intf(j));        
%         j2 = unique(j2);
%         j1 = e(j2);
%         in = ismember(j1,indelem{i});
%         elem0 = [elem0; j1(in==0)];
%         elem1 = [elem1; j1(in==1)];                                
%     end
%     elem0 = unique(elem0);
%     elem1 = unique(elem1);
    
%   find all elements that are connected to interface edgnodes
    elem0 = [];     % elements not in the process that neighbor edg nodes in the process
    elem1 = [];     % elements in the process that neighbor edg nodes NOT in the process
    for j = 1:length(intf)        
        nj = intf(j);
        rj = (mesh.cbsr_rowent2elem(nj)+1):mesh.cbsr_rowent2elem(nj+1);
        j1 = mesh.cbsr_colent2elem(rj); % list of neirghboring elements                
        in = ismember(j1,indelem{i});
        elem0 = [elem0; j1(in==0)];
        elem1 = [elem1; j1(in==1)];
    end
    elem0 = unique(elem0);
    elem1 = unique(elem1);    
    
    % [interior elements, interface elements, other elements (neighbors in other processors)]  
    elem2 = setdiff(elemNeighboringInterfEnt, [elem1; elem0]);
    if ~all(ismember(elem2,indelem{i})); error('Something wrong.'); end
    dmd{i}.elempart = [setdiff(indelem{i},[elem2; elem1]); [elem2; elem1]; elem0];
    dmd{i}.elempartpts = [length(setdiff(indelem{i},[elem2; elem1])) length(elem2)+length(elem1) length(elem0)]; 
%     elemi  = setdiff(indelem{i},elem1);
%     elemj = t2t(elem1,:);    
%     elemj = unique(elemj(:));
%     i1 = elemj==0;
%     elemj(i1) = [];
%     in = ismember(elemj,elemi);
%     elem2 = unique(elemj(in==1));    % elements in the processor that are at distance 2 to entities not in the processor     
%     dmd{i}.elempart = [setdiff(elemi,elem2); [elem2; elem1]; elem0];
%     dmd{i}.elempartpts = [length(setdiff(elemi,elem2)) length(elem2)+length(elem1) length(elem0)];        
%     dmd{i}.elempart = [setdiff(indelem{i},elem1); elem1; elem0];
%     dmd{i}.elempartpts = [length(setdiff(indelem{i},elem1)) length(elem1) length(elem0)];    

    % store elements received from neighboring subdomains to assemble linear system
    dmd{i}.elemrecv = [0*elem0 length(indelem{i})+(1:1:length(elem0))' elem0];
    for k = 1:nproc
        if k ~= i
            in = ismember(elem0,indelem{k});
            dmd{i}.elemrecv(in,1) = k;
        end
    end
    if any(dmd{i}.elemrecv(:,1) == 0); error('Domain from which element needs to be communcated was not found.'); end
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');     % This is required to sort by entries in first column       
    
    % Check that no new boundary subdomain appear for the elements compared
    % to the entities:
    newNeighSubdomains = setdiff(dmd{i}.elemrecv(:,1), dmd{i}.nbsd);
    if ~isempty(newNeighSubdomains); error('New neighboring subdomains appear for the elements compated to the entities. The DD will not work in this case'); end
    % In case this happens, it sufficies to append the new neighboring
    % subdomains to dmd{i}.nbsd (and, maybe, sort the new vector)
    
    % Check that all elements neighboring entities in processor have been
    % detected:
    elemNotDetected = setdiff(elemNeighboringEntInProc, dmd{i}.elempart);
    if ~isempty(elemNotDetected); error('Element neighboring entity in processor have not been detected.'); end
end


% store edgnodes sent to neighboring subdomains to perform the matrix-vector product
for k = 1:nproc
    dmd{k}.entsend = [];
end
for i = 1:nproc
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.entrecv(:,1)==k;
        tm = dmd{i}.entrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)
            tm(m,2) = find(dmd{k}.entpart==tm(m,3));
            if length(find(dmd{k}.entpart==tm(m,3))) ~= 1; error('Something wrong'); end
        end
        dmd{k}.entsend = [dmd{k}.entsend; tm];        
    end    
end
for k = 1:nproc
    dmd{k}.entsend = unique(dmd{k}.entsend,'rows');     % This is required to sort by entries in first column
end

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc
    dmd{k}.elemsend = [];
end
for i = 1:nproc        
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.elemrecv(:,1)==k;
        tm = dmd{i}.elemrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)
            tm(m,2) = find(dmd{k}.elempart==tm(m,3));
            if length(find(dmd{k}.elempart==tm(m,3))) ~= 1; error('Something wrong'); end
        end
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end    
end
for k = 1:nproc
    dmd{k}.elemsend = unique(dmd{k}.elemsend,'rows');     % This is required to sort by entries in first column
end


function plotpart(mesh, np, elempart, entpart)
% plot the partitioned mesh 

p = mesh.p;
t = mesh.t;

bcol = [1 1 0;... % yellow 
        1 0 1;... % magneta
        0 1 1;... % cyan
        1 0 0;... % red
        0 1 0;... % green
        0 0 1;... % blue
        1,0.4,0.6;...
        0.4,0.6,1;...
        1 1 0;... % yellow 
        1 0 1;... % magneta
        0 1 1;... % cyan
        1 0 0;... % red
        0 1 0;... % green
        0 0 1;... % blue
        1,0.4,0.6;...
        0.4,0.6,1;...
        1 1 0;... % yellow 
        1 0 1;... % magneta
        0 1 1;... % cyan
        1 0 0;... % red
        0 1 0;... % green
        0 0 1;... % blue
        1,0.4,0.6;...
        0.4,0.6,1;...
        1 1 0;... % yellow 
        1 0 1;... % magneta
        0 1 1;... % cyan
        1 0 0;... % red
        0 1 0;... % green
        0 0 1;... % blue
        1,0.4,0.6;...
        0.4,0.6,1;...
       ];
figure(2); clf;
hold on;        
for i=0:np-1
    ind = find(elempart==i);
    ti = t(ind,:);
    simpplot(p,ti,[],bcol(i+1,:));                       
end
for it=1:size(t,1)
    pmid=mean(p(t(it,:),:),1);
    txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
    text(pmid(1),pmid(2),num2str(it),txtpars{:});
end
if strcmp(mesh.hybrid,'hdg')
    for it=1:size(mesh.f,1)
        pmid=mean(p(mesh.f(it,1:2),:),1);
        i = entpart(it)+1;
        txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',bcol(7-i,:)};
        text(pmid(1),pmid(2),num2str(it),txtpars{:});
    end
elseif strcmp(mesh.hybrid,'edg')
    mesh.f(:,end+1) = 1;    
    [~,~,edg] = elconnectivities(mesh);    
    mesh.f(:,end) = [];    
    for it=1:size(edg,1)
        pmid=edg(it,:);
        i = entpart(it)+1;
        txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',bcol(7-i,:)};
        text(pmid(1),pmid(2),num2str(it),txtpars{:});
    end
end

hold off;
axis equal;      
axis tight;
axis on;  


function lia = ismemberBuiltinTypes(a,b)
numelA = numel(a);
lia = false(size(a));
b = b(:);
for i=1:numelA
    lia(i) = any(a(i)==b);   % ANY returns logical.
end
