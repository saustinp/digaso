function [extintelem,extelem] = overlappingsubdomain(intelem,t2t,overlappinglevel)

if overlappinglevel < 1
    error('overlapping level must be a positive integer');
end
    
% % global elements in nonoverlapping subdomain Omega_i
% intelem = find(elempart0==ith_sd); 

% element-to-element connectivity in nonoverlapping subdomain Omega_i
intt2t = t2t(intelem,:);

%nint = length(intelem);

% list of all elements connected to elements in nonoverlapping subdomain Omega_i
elem = unique(intt2t(:));
% remove zero or negative
if elem(1) <= 0
    elem(1) = [];
end

% exterior elements connected to elements in nonoverlapping subdomain Omega_i 
extelem = setdiff(elem,intelem);

% extent = elcon(:,extelem);
% extent = unique(extent(:));
% intfent = intersect(extent,intent);

% add more exterior elements if overlappinglevel >= 1
for n=1:overlappinglevel
    el = t2t(extelem,:);
    el = unique(el(:));   
    if el(1) <= 0
        el(1) = [];
    end
    extelem = [extelem; el];
end     
extelem = unique(extelem);
extelem = setdiff(extelem,intelem);

% while length(extelem) < length(intelem)
%     el = t2t(extelem,:);
%     el = unique(el(:));   
%     if el(1) == 0
%         el(1) = [];
%     end
%     extelem = [extelem; setdiff(el,intelem)];    
%     extelem = unique(extelem(:));
% end

% list of interior and exterior elements
extintelem = [intelem; extelem];
%extintelempts = [length(intelem) length(extelem)];

% % list of neighboring subdomains
% nbsd = unique(elempart0(extelem));




