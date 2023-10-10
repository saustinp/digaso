function [conrecv,mat] = matcon(meshp)

nproc = length(meshp);
for i = 1:nproc
    matrecv = meshp{i}.matrecv;    
    entrecv = unique(meshp{i}.cbsr_colind(matrecv+1));        
    %gentrecv = meshp{i}.entpart(entrecv);
    
    nent = length(entrecv);
    conrecv{i} = zeros(nent,nent+1);
    for j = 1:nent
        nj = entrecv(j); % local entity nj        
        rj = (meshp{i}.cbsr_rowpts(nj)+1):meshp{i}.cbsr_rowpts(nj+1);
        nbent = meshp{i}.cbsr_colind(rj); % list of its neirghboring entities            
        nbent = intersect(nbent,entrecv);
        conrecv{i}(j,2:length(nbent)+1) = [nj setdiff(nbent,nj)'];                
        %matrecv = [matrecv; meshp{i}.matrecv(j,1)*ones(length(in),1) (rj(in)-1)'];        
    end    
    nnz(conrecv{i})
    
    gentrecv = meshp{i}.entpart(entrecv);
%     for k = 1:nproc
%         if k ~= i
%             in = ismember(gentrecv,meshp{k}.intent);
%             conrecv{i}(in,1) = k;
%         end
%     end
    
    %meshp{i}.matrecv = unique(matrecv,'rows');                  
    
    n = length(matrecv);
    mat{i} = zeros(n,2);
    for j = 1:n
        r = matrecv(j);
        cn = meshp{i}.cbsr_colind(r+1);
        ii = find(meshp{i}.cbsr_rowpts<=r);        
        im = ii(end);
        jn = meshp{i}.cbsr_rowpts(im);
        rm = meshp{i}.cbsr_colind(jn+1);
        mat{i}(j,:) = [rm cn];
    end
    mat{i} = unique(mat{i},'rows');    
    
%          matrecv = [];
%         for j=1:size(meshp{i}.matrecv,1)
%             nj = meshp{i}.matrecv(j,2); % local entity nj
%             rj = (meshp{i}.cbsr_rowpts(nj)+1):meshp{i}.cbsr_rowpts(nj+1);
%             newEnt = meshp{i}.cbsr_colind(rj); % list of its neirghboring entities            
%             conEnt = intersect(newEnt,meshp{i}.matrecv(:,2));
%             [~,in] = ismember(conEnt,newEnt);
%             matrecv = [matrecv; meshp{i}.matrecv(j,1)*ones(length(in),1) (rj(in)-1)'];
%         end
%         meshp{i}.matrecv = unique(matrecv,'rows');                      
end

