function [meshp,dmd] = mkmesh_dd(mesh,nproc,preconditioner)

if nargin<=2
    preconditioner = 0;
end

elem2elem = mkt2t(mesh.t,mesh.elemtype);
if strcmp(mesh.hybrid,'hdg')
    dmd = domaindecomposition(mesh.t2f',elem2elem',nproc,preconditioner);
else
    dmd = domaindecomposition(mesh.elcon,elem2elem',nproc,preconditioner);
end
meshp = dmd;
disp('domaindecomposition complete.')

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
                % replace the global face numbering with the local face numbering
                t2f(j,k) = find(meshp{i}.entpart == t2f(j,k));
            end
        end
        for j = 1:ne
            for k = 1:nfe                
                ii = elcon(:,k,meshp{i}.elempart(j)) - (mesh.t2f(meshp{i}.elempart(j),k)-1)*npf;                
                %  on face t2f(j,k)
                ind = ((t2f(j,k)-1)*npf+1):(t2f(j,k)*npf);
                % local element-to-entity connectivities
                meshp{i}.elcon(:,k,j) = ind(ii);
            end
        end
        meshp{i}.t2f = t2f;
    else
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

%     if strcmp(mesh.hybrid,'hdg') == 0
%         matrecv = [];
%         for j=1:size(meshp{i}.matrecv,1)
% %             nj = meshp{i}.matrecv(j,3); % interface entity nj
% %             rj = (rowent2ent(nj)+1):rowent2ent(nj+1);
% %             newEnt = colent2ent(rj); % list of its neirghboring entities            
% %             conEnt = intersect(newEnt,meshp{i}.matrecv(:,3));
% %             [~,in] = ismember(conEnt,meshp{i}.matrecv(:,3));
% %             locEnt = unique(meshp{i}.matrecv([j;in],2));            
%             nj = meshp{i}.matrecv(j,2); % local entity nj
%             rj = (meshp{i}.cbsr_rowpts(nj)+1):meshp{i}.cbsr_rowpts(nj+1);
%             newEnt = meshp{i}.cbsr_colind(rj); % list of its neirghboring entities            
%             conEnt = intersect(newEnt,meshp{i}.matrecv(:,2));
%             [~,in] = ismember(conEnt,newEnt);            
% %             newEnt
% %             conEnt            
% %             in
% %             nj
% %             meshp{i}.cbsr_colind(rj(in))
% %             rj(in)
% %             pause            
%             matrecv = [matrecv; meshp{i}.matrecv(j,1)*ones(length(in),1) (rj(in)-1)'];
%         end
%         meshp{i}.matrecv = unique(matrecv,'rows');
%         %pause
%     end
%     
%     for j = 1:length(meshp{i}.nbsd)
%         meshp{i}.matsendpts(j) = length(find(meshp{i}.matsend(:,1)==meshp{i}.nbsd(j)));
%         meshp{i}.matrecvpts(j) = length(find(meshp{i}.matrecv(:,1)==meshp{i}.nbsd(j)));
%     end    
%     meshp{i}.matsend(:,2) = meshp{i}.cbsr_rowpts(meshp{i}.matsend(:,2));
%     %meshp{i}.matrecv(:,2) = meshp{i}.cbsr_rowpts(meshp{i}.matrecv(:,2));   
%     meshp{i}.matsend = meshp{i}.matsend(:,2);
%     meshp{i}.matrecv = meshp{i}.matrecv(:,2);        
end
disp('Postprocessing No. 1 complete.')

if strcmp(mesh.hybrid,'hdg') 
    for i = 1:nproc
        for j = 1:length(meshp{i}.nbsd)
            meshp{i}.matsendpts(j) = length(find(meshp{i}.matsend(:,1)==meshp{i}.nbsd(j)));
            meshp{i}.matrecvpts(j) = length(find(meshp{i}.matrecv(:,1)==meshp{i}.nbsd(j)));
        end    
        meshp{i}.matsend(:,2) = meshp{i}.cbsr_rowpts(meshp{i}.matsend(:,2));
        meshp{i}.matrecv(:,2) = meshp{i}.cbsr_rowpts(meshp{i}.matrecv(:,2));   
        meshp{i}.matsend = meshp{i}.matsend(:,2);
        meshp{i}.matrecv = meshp{i}.matrecv(:,2);       
    end
else
    % store faces sent to neigboring cpus to assemble the preconditioner
    for i = 1:nproc
         meshp{i}.matsend = [];
         matrecv = [];
        for j=1:size(meshp{i}.matrecv,1)
            nj = meshp{i}.matrecv(j,2); % local entity nj
            rj = (meshp{i}.cbsr_rowpts(nj)+1):meshp{i}.cbsr_rowpts(nj+1);
            newEnt = meshp{i}.cbsr_colind(rj); % list of its neighboring entities            
            % only take neighbors in meshp{i}.matrecv(:,2)
            conEnt = intersect(newEnt,meshp{i}.matrecv(:,2)); 

%         nj = dmd{i}.matrecv(j,3); % exterior entity nj
%         rj = (rowent2elem(nj)+1):rowent2elem(nj+1);
%         ei = colent2elem(rj); % list of its neighboring elements                
            
            [~,in] = ismember(conEnt,newEnt);
            matrecv = [matrecv; meshp{i}.matrecv(j,1)*ones(length(in),1) (rj(in)-1)'];
        end
        meshp{i}.matrecv = unique(matrecv,'rows');                  
    end
    for i = 1:nproc          
        for j = 1:length(meshp{i}.nbsd)
            % cpu k sends information to cpu i
            k = meshp{i}.nbsd(j);
            ii = meshp{i}.matrecv(:,1)==k;
            tm = meshp{i}.matrecv(ii,:);
            tm(:,1) = i;
            meshp{k}.matsend = [meshp{k}.matsend; tm];        
        end    
    end
    for i = 1:nproc   
        meshp{i}.matsend = unique(meshp{i}.matsend,'rows');            
        for j = 1:length(meshp{i}.nbsd)
            meshp{i}.matrecvpts(j) = length(find(meshp{i}.matrecv(:,1)==meshp{i}.nbsd(j)));
            meshp{i}.matsendpts(j) = length(find(meshp{i}.matsend(:,1)==meshp{i}.nbsd(j)));        
        end        
        meshp{i}.matrecv = meshp{i}.matrecv(:,2);        
        meshp{i}.matsend = meshp{i}.matsend(:,2);
    end
end
disp('Postprocessing No. 2 complete.')
