function dmd = domaindecompositionmpi(elcon,t2t,t2f,nproc,preconditioner,elconhdg)
% elcon: element-to-entity connectivity
% t2t: element-to-element connectivity  
% t2f: element-to-face connectivity  
% nproc: number of processors
% preconditioner: 0 for block Jacobi and 1-entity RAS
%                 1 for 1-element RAS
%                 2 for 2-element RAS
% elconhdg: element-to-entity connectivity for hdg
% dmd: domain decomposition structure

% determine if this is hdg
if max(abs(size(t2f')-size(elcon)))==0
    hybrid = 'hdg';
    nfe = size(t2f,2);
    npf = size(elconhdg,1)/nfe;
    elconhdg = reshape(elconhdg,[npf nfe size(elcon,2)]);
else
    hybrid = 'edg';
end

if nargin<=3
    preconditioner = 0;
end
overlappinglevel = preconditioner+2; % determine the level of overlap between subdomains

ne = size(elcon,2);   % number of elements
ndof = max(elcon(:)); % number of entities

if max(abs((1:ndof)'-unique(elcon(:))))
    error('It looks like there are ghost entities in the mesh.'); 
end

% nonoverlapping domain decomposition using METIS
[elemall, entall] = elconpartition(elcon,nproc);
elemall = elemall+1;
entall = entall+1;

% loop over each subdomain
checkElements = 0;
checkEntities = 0;
for i=1:nproc
    % list of  elements in nonoverlapping subdomain i                       
    dmd{i}.intelem = find(elemall==i); 
    checkElements = checkElements + length(dmd{i}.intelem);
    
    % list of  entities in nonoverlapping subdomain i                       
    dmd{i}.intent = find(entall==i); 
    checkEntities = checkEntities + length(dmd{i}.intent);
end

numElem = ne;
numEntities = ndof;
if checkElements ~= numElem; error('Domain decomposition for elements was not performed properly.'); end
if checkEntities ~= numEntities; error('Domain decomposition for entities was not performed properly.'); end

% Check that all elements in the processor have at least one entity in the
% processor
disp(' ');
disp('Check that all elements in the processor are connected to at least one entity in the processor...');
for i=1:nproc
    disp(['Proc No. ', num2str(i), ' / ', num2str(nproc)]);
    elementsToRemove = [];
    for j=1:length(dmd{i}.intelem)
        if isempty(intersect(dmd{i}.intent,elcon(:,dmd{i}.intelem(j)))) 
            entitiesInElement = elcon(:,dmd{i}.intelem(j));
            candidateProcessors = entall(entitiesInElement);
            k = mode(candidateProcessors);
            dmd{k}.intelem = sort([dmd{k}.intelem;dmd{i}.intelem(j)]);
            elementsToRemove = [elementsToRemove;dmd{i}.intelem(j)];
            warning('An element in the processor did not have any entity in the processors. The issue has been fixed successfully.');
        end
    end
    if ~isempty(elementsToRemove)
        dmd{i}.intelem = sort(setdiff(dmd{i}.intelem,elementsToRemove));
    end
end

disp(' ');
disp('Check that all elements in the processor are connected to at least one element in the processor...');
for i=1:nproc
    disp(['Proc No. ', num2str(i), ' / ', num2str(nproc)]);
    elementsInProcessor = t2t(dmd{i}.intelem,:);
    if ~isempty(setdiff(dmd{i}.intelem,unique(elementsInProcessor(:))))
        error('Some elements in the processor are not connected to any element in the processor.');
    end
end

disp(' ');
disp('Check that all entities in the processor are connected to at least one element in the processor ...');
for i=1:nproc
    disp(['Proc No. ', num2str(i), ' / ', num2str(nproc)]);
    entitiesInProcessor = elcon(:,dmd{i}.intelem);
    if ~isempty(setdiff(dmd{i}.intent,unique(entitiesInProcessor(:))))
        error('Some entities in the processor are not connected to any element in the processor.');
    end
end
disp(' ');

% (rowent2elem,colent2elem) entity-to-element connectivity
% (rowent2ent,colent2ent) entity-to-entity connectivity
%[rowent2elem,colent2elem,rowent2ent,colent2ent] = elcon2entcon(elcon);

% loop over each subdomain
for i=1:nproc    
    disp(['Preprocessing subdomain No. ',num2str(i)]);
    
    % list of elements on overlapping subdomain i
    [extintelem,~] = overlappingsubdomain(dmd{i}.intelem,t2t,overlappinglevel);        
        
    % (rowent2elem,colent2elem) entity-to-element connectivity
    % (rowent2ent,colent2ent) entity-to-entity connectivity        
    [rowent2elem,colent2elem,rowent2ent,colent2ent,~,ent2ind] = elcon2entconmpi(elcon(:,extintelem),extintelem);
    
    dmd{i}.rowent2elem = rowent2elem;
    dmd{i}.colent2elem = colent2elem;
    dmd{i}.rowent2ent  = rowent2ent;
    dmd{i}.colent2ent = colent2ent;
    dmd{i}.ent2ind = ent2ind;
    
    % Find all elements connected to entities in processor (i.e. elements
    % required to compute the rows of entities in the processor)
    elemNeighboringEntInProc = [];
    for j = 1:length(dmd{i}.intent)
        nj = dmd{i}.intent(j);        
        ij = ent2ind(nj);
        rj = (rowent2elem(ij)+1):rowent2elem(ij+1);
        newElem = colent2elem(rj);
        elemNeighboringEntInProc = [elemNeighboringEntInProc; newElem(:)]; 
    end
    elemNeighboringEntInProc = unique(elemNeighboringEntInProc(:));
    
    % find all interface entities in dmd{i}.intent.
    nent = length(dmd{i}.intent);
    entnumber = zeros(length(dmd{i}.intent),1);    
    for j = 1:nent  % For each entity in in dmd{i}.intent
        nj = dmd{i}.intent(j);  % entity nj
        ij = ent2ind(nj);
        rj = (rowent2ent(ij)+1):rowent2ent(ij+1); % positions of row nj 
        entnb = colent2ent(rj); % list of its neirghboring entities 
        if all(ismember(entnb,dmd{i}.intent))            
            entnumber(j) = 2; % entity nj is inside the subdomain i 
        else            
            entnumber(j) = 1; % entity nj is on the interface between two subdomains            
        end                
    end
    
    % list of all interface edgnodes in the subdomain
    intfent = dmd{i}.intent(entnumber==1);

    % Find all elements connected to interface entities (i.e. elements
    % required to compute the rows of interface)
    elemNeighboringInterfEnt = [];
    for j = 1:length(intfent)
        nj = intfent(j); % interface entity nj
        ij = ent2ind(nj);
        rj = (rowent2elem(ij)+1):rowent2elem(ij+1);
        newElem = colent2elem(rj); % list of its neirghboring elements
        elemNeighboringInterfEnt = [elemNeighboringInterfEnt; newElem(:)]; 
    end
    elemNeighboringInterfEnt = unique(elemNeighboringInterfEnt(:));
    
    % Find all entities connected to interface entities
    entNeighboringInterfEnt = [];
    for j = 1:length(intfent)
        nj = intfent(j); % interface entity nj
        ij = ent2ind(nj);
        rj = (rowent2ent(ij)+1):rowent2ent(ij+1);
        newEnt = colent2ent(rj); % list of its neirghboring entities
        entNeighboringInterfEnt = [entNeighboringInterfEnt; newEnt(:)]; 
    end
    entNeighboringInterfEnt = unique(entNeighboringInterfEnt(:));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % interface elements
    intfelem = intersect(elemNeighboringInterfEnt,dmd{i}.intelem);
    % exterior elements
    extelem = setdiff(elemNeighboringInterfEnt,intfelem);
    % n-element RAS preconditioner
    if preconditioner>0
        % 1-element RAS preconditioner
        el = t2t(intfelem,:); % find neighboring elements of the interface elements
        el = unique(el(:)); % remove duplications  
        el = el(el > 0);    % remive zero indices
        % add only exterior neighboring elements
        extelem = [extelem; setdiff(dmd{i}.intelem,el)];
        % repeat the process for n-element RAS preconditioner
        for n=2:1:preconditioner
            el = t2t(extelem,:);
            el = unique(el(:));   
            el = el(el > 0);
            extelem = [extelem; setdiff(dmd{i}.intelem,el)];
        end         
    end
    extelem = unique(extelem(:));
        
    % [interior elements, interface elements, exterior elements]  
    dmd{i}.elempart = [setdiff(dmd{i}.intelem,intfelem); intfelem; extelem];
    dmd{i}.elempartpts = [length(dmd{i}.elempart)-(length(intfelem)+length(extelem)) length(intfelem) length(extelem)]; 
    % interior elements are elements which are in intelem and not connected to interface entities  
    % interface elements are elements which are in intelem and connected to interface entities  
    % exterior elements are elements which are not in intelem. 
    % exterior elements depends on the preconditioner parameter. At the
    % very least they should include all the elements that are connected to
    % interface entities, so that elemNeighboringEntInPro is a subset of elempart. 
    %  - For Block-Jacobi and 1-entity RAS: it is the least requirement.
    %  - For 1-element RAS: it is the least requirement plus the neighboring elements
    %                       of interface elements
    %  - For 2-element RAS: it is the 1-element RAS plus the neighboring
    %                       elements of the 1-element RAS
    
%     dmd{i}.elempart'
%     dmd{i}.elempartpts
    
    % store elements received from neighboring subdomains to assemble linear system
    dmd{i}.elemrecv = [0*extelem length(dmd{i}.intelem)+(1:1:length(extelem))' extelem];
    for k = 1:nproc
        if k ~= i
            in = ismember(extelem,dmd{k}.intelem);
            dmd{i}.elemrecv(in,1) = k;
        end
    end    
    if any(dmd{i}.elemrecv(:,1) == 0) 
        error('Domain from which element needs to be communcated was not found.');
    end
     % This is required to sort by entries in first column           
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');    
    %dmd{i}.elemrecv
    
    % find all entities that are connected to interface elements and exterior elements
    en = elcon(:,[intfelem; extelem]);
    en = unique(en(:));
    en = en(en>0);
    % exterior entities
    extent = setdiff(en,dmd{i}.intent);
            
    %extent'
    %length(extent)
    
    % [interior entities, interface entities, exterior entities]
    dmd{i}.entpart = [setdiff(dmd{i}.intent,intfent); intfent; extent];
    dmd{i}.entpartpts = [length(dmd{i}.entpart)-(length(intfent)+length(extent)) length(intfent) length(extent)];        
    % interior entities are entities which are in intent and not connected to interface entities  
    % interface entities are entities which are in intent and have exterior entities as neighbors.
    % exterior entities are entities which are not in intent. 
    % exterior entities are are connected to interface elements and exterior elements.
    % hence, exterior entities depends on the preconditioner parameter. 
    
    %dmd{i}.entpart'
    %dmd{i}.entpartpts
    
    % find and store exterior entities received from neighboring subdomains to 
    % apply the preconditioner to a vector, because the self processor does not 
    % compute the block vectors associated with the exterior entities. They have
    % to be computed by neighboring processors and sent to the self processor. 
    % [neighboring subdomains, local numberings, global numberings]
    dmd{i}.entrecv = [0*extent length(dmd{i}.intent)+(1:1:length(extent))' extent];    
    for k = 1:nproc
        if k ~= i
            in = ismember(extent,dmd{k}.intent);
            dmd{i}.entrecv(in,1) = k;
        end
    end
    if any(dmd{i}.entrecv(:,1) == 0)
        error('Domain from which entity needs to be communcated was not found.'); 
    end
    % It is required to sort by entries in the first column  
    dmd{i}.entrecv = unique(dmd{i}.entrecv,'rows');       
    % list of neighboring subdomains from which info is required for matrix-vector product    
    dmd{i}.nbsd = unique(dmd{i}.entrecv(:,1))';             
        
    %dmd{i}.entrecv
    %dmd{i}.nbsd
    
    % find exterior entities connected to interface entities to perform matrix-vector product. 
    on = [];
    for j = 1:size(dmd{i}.entrecv,1)
        nj = dmd{i}.entrecv(j,3); % exterior entity nj
        ij = ent2ind(nj);        
        rj = (rowent2ent(ij)+1):rowent2ent(ij+1);
        ei = colent2ent(rj); % list of its neighboring entities
        if any(ismember(ei,intfent)) 
            % any neighboring entities of the entity nj are inside the interface entities
            % entity nj is connected to interface entities
            on = [on j];                    
        end
    end            
    dmd{i}.vecrecv = dmd{i}.entrecv(on,:);
    % It is required to sort by entries in the first column  
    dmd{i}.vecrecv = unique(dmd{i}.vecrecv,'rows');       
    
    %dmd{i}.vecrecv
    
    % find exterior entities for which any of neighboring elements is not in elempart.
    % Because the block matrices associated with those entities can not be
    % formed by the self processor. They have to be computed by neighboring
    % processors and sent to the self processor. matrecv is only needed for RAS, 
    % since BJ require no communication.     
    % find and store entities received from neigboring subdomains to assemble the preconditioner    
    %dmd{i}.matrecv = dmd{i}.entrecv;
    on = [];
    for j = 1:size(dmd{i}.entrecv,1)
        nj = dmd{i}.entrecv(j,3); % exterior entity nj
        ij = ent2ind(nj);
        rj = (rowent2elem(ij)+1):rowent2elem(ij+1);
        ei = colent2elem(rj); % list of its neighboring elements                        
        %if all(ismember(ei,dmd{i}.elempart)) 
        if isempty(setdiff(ei,dmd{i}.elempart))==0   
            % one of the neighboring elements of the entity nj is NOT inside the overlapped subdomain i             
            on = [on j];                    
        end
    end            
    dmd{i}.matrecv = dmd{i}.entrecv(on,:);
    % It is required to sort by entries in the first column  
    dmd{i}.matrecv = unique(dmd{i}.matrecv,'rows');       
    
%     on
%     dmd{i}.matrecv
%     size(dmd{i}.matrecv)
    
    if strcmp(hybrid,'hdg')        
        dmd{i}.elcon = elconhdg(:,:,dmd{i}.elempart);
        ne = length(dmd{i}.elempart);
        t2fhdg = t2f(dmd{i}.elempart,:);
        for j = 1:ne
            for k = 1:nfe
                % replace the global face numbering with the local face numbering
                t2fhdg(j,k) = find(dmd{i}.entpart == t2fhdg(j,k));
            end
        end
        for j = 1:ne
            for k = 1:nfe                
                ii = elconhdg(:,k,dmd{i}.elempart(j)) - (t2f(dmd{i}.elempart(j),k)-1)*npf;                
                %  on face t2f(j,k)
                ind = ((t2fhdg(j,k)-1)*npf+1):(t2fhdg(j,k)*npf);
                % local element-to-entity connectivities
                dmd{i}.elcon(:,k,j) = ind(ii);
            end
        end
        dmd{i}.t2f = t2fhdg;
        dmd{i}.elcon = reshape(dmd{i}.elcon,[npf*nfe ne]);   
        
        % local entity-to-element and entity-to-entity connectivities
        [dmd{i}.bcrs_rowent2elem,dmd{i}.bcrs_colent2elem,...
         dmd{i}.bcrs_rowent2ent,dmd{i}.bcrs_colent2ent] = elcon2entcon(dmd{i}.t2f');            
    else
        % global element-to-entity connectivity
        dmd{i}.elcon = elcon(:,dmd{i}.elempart); 
        % mapping from global to local
        entMapping = zeros(max(dmd{i}.entpart),1);
        entMapping(dmd{i}.entpart) = 1:length(dmd{i}.entpart);
        % local element-to-entity connectivity
        dmd{i}.elcon = entMapping(dmd{i}.elcon);
        if min(dmd{i}.elcon(:)) < 1; error('Something wrong.'); end
        
        %dmd{i}.elcon
        
        % global element-to-face connectivity
        dmd{i}.t2f = t2f(dmd{i}.elempart,:);
        % global faces in the subdomain
        facesInProcessor = dmd{i}.t2f(:);    
        facesInProcessor = unique(facesInProcessor);
        facesInProcessor = facesInProcessor(facesInProcessor > 0);
                
        %dmd{i}.t2f
        %facesInProcessor'                
        
        % mapping from global to local
        faceMapping = zeros(max(facesInProcessor),1);
        faceMapping(facesInProcessor) = 1:length(facesInProcessor);
        % local element-to-face connectivity
        dmd{i}.t2f = faceMapping(dmd{i}.t2f);       
        if min(dmd{i}.t2f(:)) < 1; error('Something wrong.'); end
                
        % local entity-to-element and entity-to-entity connectivities
        [dmd{i}.bcrs_rowent2elem,dmd{i}.bcrs_colent2elem,...
         dmd{i}.bcrs_rowent2ent,dmd{i}.bcrs_colent2ent] = elcon2entcon(dmd{i}.elcon);                    
    end
    
    tm = dmd{i}.bcrs_rowent2ent(2:end)-dmd{i}.bcrs_rowent2ent(1:end-1);
    dmd{i}.maxBlocksPerRow = max(tm);
    dmd{i}.minBlocksPerRow = min(tm);
    dmd{i}.bcrs_nrows = length(dmd{i}.bcrs_rowent2ent)-1;
    dmd{i}.bcrs_nblks = length(dmd{i}.bcrs_colent2ent);    
    
    numTotalEntities = max(dmd{i}.elcon(:));
    numRealEntities = length(unique(dmd{i}.elcon(:)));    
    if numRealEntities ~= numTotalEntities
        error('It looks like there are ghost (missing) entities in the mesh.');       % The construction of BJ and BK from cbsr is based on the assumption that there are no ghost entities
    end
    
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
    
    entNotDetected = setdiff(entNeighboringInterfEnt, dmd{i}.entpart);
    if ~isempty(entNotDetected); error('Entity neighboring entity in processor have not been detected.'); end        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
end

% store entities sent to neighboring subdomains to perform the matrix-vector product
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
  % This is required to sort by entries in first column
    dmd{k}.entsend = unique(dmd{k}.entsend,'rows');   
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
 % This is required to sort by entries in first column
 for k = 1:nproc
    dmd{k}.elemsend = unique(dmd{k}.elemsend,'rows');    
 end

 % store entities sent to neighboring subdomains to perform the matrix-vector product
for k = 1:nproc
    dmd{k}.vecsend = [];
end
for i = 1:nproc
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j); 
        ii = dmd{i}.vecrecv(:,1)==k;
        tm = dmd{i}.vecrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)
            tm(m,2) = find(dmd{k}.entpart==tm(m,3));
            if length(find(dmd{k}.entpart==tm(m,3))) ~= 1; error('Something wrong'); end
        end
        dmd{k}.vecsend = [dmd{k}.vecsend; tm];        
    end    
end
for k = 1:nproc
  % This is required to sort by entries in first column
    dmd{k}.vecsend = unique(dmd{k}.vecsend,'rows');   
end

% store faces sent to neigboring cpus to assemble the preconditioner
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
% This is required to sort by entries in first column
for k = 1:nproc   
    dmd{k}.matsend = unique(dmd{k}.matsend,'rows');        
end

if strcmp(hybrid,'hdg') 
    for i = 1:nproc
        dmd{i}.matsend(:,2) = dmd{i}.bcrs_rowent2ent(dmd{i}.matsend(:,2));
        dmd{i}.matrecv(:,2) = dmd{i}.bcrs_rowent2ent(dmd{i}.matrecv(:,2));           
    end
else    
    irecv = ones(nproc,1);
    isend = ones(nproc,1);
    for i = 1:nproc
        matrecv{i} = dmd{i}.matrecv;            
        matsend{i} = dmd{i}.matsend;            
        nent = size(matrecv{i},1);        
        dmd{i}.matrecv = zeros(nent*dmd{i}.maxBlocksPerRow,2);                
        nent = size(matsend{i},1);        
        dmd{i}.matsend = zeros(nent*dmd{i}.maxBlocksPerRow,2);                
    end
    for i = 1:nproc
        
        % list of elements on overlapping subdomain i
        [extintelem,~] = overlappingsubdomain(dmd{i}.intelem,t2t,overlappinglevel);
        
        % (rowent2elem,colent2elem) entity-to-element connectivity
        % (rowent2ent,colent2ent) entity-to-entity connectivity    
        [rowent2elem,colent2elem,~,~,~,ent2ind] = elcon2entconmpi(elcon(:,extintelem),extintelem);
        
        nent = size(matrecv{i},1);                
        for j = 1:nent
            ncpu = matrecv{i}(j,1); % neighboring cpu
            nj = matrecv{i}(j,2);   % local entity nj          
            rj = (dmd{i}.bcrs_rowent2ent(nj)+1):dmd{i}.bcrs_rowent2ent(nj+1);
            % list of local neighboring entities to nj            
            nbent = dmd{i}.bcrs_colent2ent(rj); 
            % take neighbors in matrecv{i}(:,2)
            nbentj = intersect(nbent,matrecv{i}(:,2));
            % reorder the list of neighboring entities
            nbentj = [nj setdiff(nbentj,nj)'];      
            [~,ij] = ismember(nbentj,nbent);
            % get a subset of rj
            rowj = rj(ij);            
%             if j==1
%                 [i-1 ij]
%                 [i-1 nbent']
%                 [i-1 nbentj]                
%                 [i-1 rowj]
%             end
            
            % add [ncpu rowj(1)-1] to matrecv
            dmd{i}.matrecv(irecv(i),:) = [ncpu rowj(1)-1];
            irecv(i) = irecv(i) + 1;                        
            
            gj = dmd{i}.entpart(nj); % global entity gj 
            ij = ent2ind(gj);
            rj = (rowent2elem(ij)+1):rowent2elem(ij+1);
            % list of global neighboring elements to gj
            nbelemj = colent2elem(rj);
            
%             if j==1           
%                 %[i-1 ncpu gj ij length(rj)]
%                 [i-1 nbelemj']                                
%             end            
            
            % find an element in matsend{ncpu}(:,3) to match gj
            m = matsend{ncpu}(:,3) == gj; 
            % local entity mj on neighboring ncpu              
            mj = matsend{ncpu}(m,2); 
            pj = (dmd{ncpu}.bcrs_rowent2ent(mj(1))+1):dmd{ncpu}.bcrs_rowent2ent(mj(1)+1);
            % add [i pj(1)-1] to matsend on neighboring ncpu                      
            dmd{ncpu}.matsend(isend(ncpu),:) = [i pj(1)-1];
            isend(ncpu) = isend(ncpu)+1;                                                    
            % list of local neighboring entities to mj            
            mbent = dmd{ncpu}.bcrs_colent2ent(pj); 
            % list of global neighboring entities to gj on neighboring ncpu                      
            gbent = dmd{ncpu}.entpart(mbent); 
            
%             if j==1                                                             
%                 gbent'
%             end
            
            for k = 2:length(nbentj)
                gk = dmd{i}.entpart(nbentj(k)); % global entity gk 
                ik = ent2ind(gk);
                rk = (rowent2elem(ik)+1):rowent2elem(ik+1);
                
                %[gk ik length(rk) rowent2elem(ik) rowent2elem(ik+1)]
                
                % list of global neighboring elements to gk
                nbelemk = colent2elem(rk);                                
                
                % intersect nbelemj and nbelemk to get a common set
                nbelem = intersect(nbelemj,nbelemk);
                                        
                % if there are more than one common element
                % if length(nbelem)>1 
                % if nbelem is NOT a subset of elempart
                if isempty(setdiff(nbelem,dmd{i}.elempart))==0
                    % add [ncpu rowj(k)-1] to matrecv
                    dmd{i}.matrecv(irecv(i),:) = [ncpu rowj(k)-1];
                    irecv(i) = irecv(i) + 1; 
                    % find an element in gbent to match gk
                    nk = gbent==gk;
                    % add [i pj(nk)-1] to matsend
                    dmd{ncpu}.matsend(isend(ncpu),:) = [i pj(nk)-1];
                    isend(ncpu) = isend(ncpu)+1;                    
                end
            end            
        end    
    end    
    for i = 1:nproc
        dmd{i}.matrecv = dmd{i}.matrecv(1:(irecv(i)-1),:);
        dmd{i}.matsend = dmd{i}.matsend(1:(isend(i)-1),:);        
    end    
end

for i = 1:nproc    
    for j = 1:length(dmd{i}.nbsd)
        dmd{i}.entsendpts(j) = length(find(dmd{i}.entsend(:,1)==dmd{i}.nbsd(j)));
        dmd{i}.entrecvpts(j) = length(find(dmd{i}.entrecv(:,1)==dmd{i}.nbsd(j)));        
        dmd{i}.elemsendpts(j) = length(find(dmd{i}.elemsend(:,1)==dmd{i}.nbsd(j)));
        dmd{i}.elemrecvpts(j) = length(find(dmd{i}.elemrecv(:,1)==dmd{i}.nbsd(j)));
        dmd{i}.vecsendpts(j) = length(find(dmd{i}.vecsend(:,1)==dmd{i}.nbsd(j)));
        dmd{i}.vecrecvpts(j) = length(find(dmd{i}.vecrecv(:,1)==dmd{i}.nbsd(j)));                
        dmd{i}.matsendpts(j) = length(find(dmd{i}.matsend(:,1)==dmd{i}.nbsd(j)));
        dmd{i}.matrecvpts(j) = length(find(dmd{i}.matrecv(:,1)==dmd{i}.nbsd(j)));        
    end
    dmd{i}.entsend = dmd{i}.entsend(:,2);
    dmd{i}.entrecv = dmd{i}.entrecv(:,2);
    dmd{i}.elemsend = dmd{i}.elemsend(:,2);
    dmd{i}.elemrecv = dmd{i}.elemrecv(:,2);
    dmd{i}.vecsend = dmd{i}.vecsend(:,2);
    dmd{i}.vecrecv = dmd{i}.vecrecv(:,2);
    dmd{i}.matsend = dmd{i}.matsend(:,2);
    dmd{i}.matrecv = dmd{i}.matrecv(:,2);           
end

