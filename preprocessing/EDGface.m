
function isEDGface = EDGface(p, t, f, nhdgfaces, requiredHDGfaces, auxInputs, methodFlag)

% HDG faces are the union of those given by the method specified by "methodFlag"
% and those in the array requiredHDGfaces.
% methodFlag = 1: Decide faces based on minimum face heights.
% methodFlag = 2: Decide faces based on (faces neighboring) minimum element
%                 height. auxInputs = t2f.
% methodFlag = 3: Decide faces based on (faces neighboring) given list of element
%                 sizes. auxInputs = {t2f, he}.

if nargin < 4; nhdgfaces = 0.1; end % HDG faces take 10 percent 
if nargin < 5; requiredHDGfaces = []; end
if nargin < 7; methodFlag = 1; end

nf = size(f,2);

if methodFlag == 1       % Decide HDG faces based on face size
    hf = facesize(p,f);
    hf = min(hf,[],1); % face sizes
     
    % sort faces according to their sizes
    [~,ind] = sort(hf);
    ind = ind(:);

    if nhdgfaces < 1
        % number of hdgfaces
        nhdgfaces = round(nhdgfaces*nf);
    end
    nhdgfaces = min(nf, nhdgfaces);
    ihdg = ind(1:nhdgfaces); % indices for HDG faces

elseif methodFlag == 2               % Decide HDG faces based on element size
    if nargin < 6; error('t2f needs to be provided as input argument.'); end
    t2f = auxInputs;
    
    nhdgelements = nhdgfaces;
    
    he = meshsize(p,t);
    he = min(he,[],1); % element sizes
    
    [~,ind] = sort(he);
    ind = ind(:);
    
    ne = size(t,2);
    if nhdgelements < 1    
        nhdgelements = round(nhdgelements*ne);
    end
    nhdgelements = min(ne, nhdgelements);
    ehdg = ind(1:nhdgelements); % indices for HDG elements
    
    ihdg = t2f(ehdg,:);
    ihdg = unique(ihdg(:)); % indices for HDG faces
elseif methodFlag == 3               % Decide HDG faces based on given element sizes
    if nargin < 6; error('{t2f,he} need to be provided as input argument.'); end
    if length(auxInputs) ~= 2; error('{t2f,he} needs to be provided as input argument.'); end
    
    t2f = auxInputs{1};
    he = auxInputs{2};
    
    nhdgelements = nhdgfaces;
    
    [~,ind] = sort(he);
    ind = ind(:);
    
    ne = size(t,2);
    if nhdgelements < 1    
        nhdgelements = round(nhdgelements*ne);
    end
    nhdgelements = min(ne, nhdgelements);
    ehdg = ind(1:nhdgelements); % indices for HDG elements
    
    ihdg = t2f(ehdg,:);
    ihdg = unique(ihdg(:)); % indices for HDG faces
else
    error('Flag type not recognized.')
end

% Include required HDG faces:
if ~isempty(requiredHDGfaces)
    ihdg = [ihdg(:); requiredHDGfaces(:)]; % indices for HDG faces
    ihdg = sort(ihdg);
    ihdg = unique(ihdg(:));
end

isEDGface = ones(nf,1); % EDG faces
isEDGface(ihdg) = 0; % HDG faces 

plotFlag = 0;
if plotFlag == 1
    pHDG = f(:,ihdg);
    pEDG = setdiff(1:size(p,1),pHDG(:));
    figure();
    plot(p(pHDG(:),1),p(pHDG(:),2),'.');
    hold on
    plot(p(pEDG(:),1),p(pEDG(:),2),'.');
end

end


% function isEDGface = EDGface(p,t,f,nhdgfaces)
% 
% if nargin<4
%     nhdgfaces = 0.1; % HDG faces take 10 percent 
% end
% 
% nd = size(p,2);
% [nvf,nf] = size(f);
% 
% pf = reshape(p(f,:),[nvf nf nd]);
% if nvf == 2
%     v1 = 1;
%     v2 = 2;
% elseif nvf == 3
%     v1 = [1 2 3];
%     v2 = [2 3 1];    
% elseif nvf == 4
%     v1 = [1 2 3 4];
%     v2 = [2 3 4 1];    
% end
% 
% nv = length(v1);
% hf = zeros(nv,nf);
% for j = 1:nv
%     i1 = v1(j);
%     i2 = v2(j);
%     hf(j,:) = (pf(i1,:,1)-pf(i2,:,1)).^2;
%     for k = 2:nd
%         hf(j,:) = hf(j,:) + (pf(i1,:,k)-pf(i2,:,k)).^2;
%     end    
% end
% hf = sqrt(hf); 
% hf = min(hf,[],1); % face sizes
% 
% % sort faces according to their sizes
% [~,ind] = sort(hf);
% 
% if nhdgfaces < 1
%     % number of hdgfaces
%     nhdgfaces = round(nhdgfaces*nf);
% end
% ihdg = ind(1:nhdgfaces); % indices for HDG faces
% 
% isEDGface = ones(nf,1); % EDG faces
% isEDGface(ihdg) = 0; % HDG faces 

function hf = facesize(p,f)

nd = size(p,2);
[nvf,nf] = size(f);

pf = reshape(p(f,:),[nvf nf nd]);
if nvf == 2
    v1 = 1;
    v2 = 2;
elseif nvf == 3
    v1 = [1 2 3];
    v2 = [2 3 1];    
elseif nvf == 4
    v1 = [1 2 3 4];
    v2 = [2 3 4 1];    
end

nv = length(v1);
hf = zeros(nv,nf);
for j = 1:nv
    i1 = v1(j);
    i2 = v2(j);
    hf(j,:) = (pf(i1,:,1)-pf(i2,:,1)).^2;
    for k = 2:nd
        hf(j,:) = hf(j,:) + (pf(i1,:,k)-pf(i2,:,k)).^2;
    end    
end
hf = sqrt(hf); 
hf = min(hf,[],1); % face sizes

end


function he = meshsize(p,t)

nd = size(p,2);
[nve,ne] = size(t);
pe = reshape(p(t,:),[nve ne nd]);

if nve == 2
    v1 = 1;
    v2 = 2;
elseif nve == 3
    v1 = [1 2 3];
    v2 = [2 3 1];    
elseif nve == 4 && nd==2  % quad element
    v1 = [1 2 3 4];
    v2 = [2 3 4 1];    
elseif nve == 4 && nd == 3  % tet element
    v1 = [1 2 3 4 4 4];
    v2 = [2 3 1 1 2 3];        
elseif nve == 8
    v1 = [1 2 3 4 5 6 7 8 1 2 3 4];
    v2 = [2 3 4 1 6 7 8 5 5 6 7 8];
else
    error('Element type not recognized')
end

nv = length(v1);
he = zeros(nv,ne);
for j = 1:nv
    i1 = v1(j);
    i2 = v2(j);
    he(j,:) = (pe(i1,:,1)-pe(i2,:,1)).^2;
    for k = 2:nd
        he(j,:) = he(j,:) + (pe(i1,:,k)-pe(i2,:,k)).^2;
    end    
end
he = sqrt(he);

end
