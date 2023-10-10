function ent2cpu = edgpartition(facecon, t, t2f, re, ce, elem2cpu, face2cpu, nproc, overlappinglevel)

entmax = max(facecon(:))+1;
ent2cpu = -ones(entmax,1);            
nent = 0;
for i = 1:nproc
    % list of  elements in nonoverlapping subdomain i                       
    intelem = find(elem2cpu==(i-1));     

    % list of  faces in nonoverlapping subdomain i                       
    intface = find(face2cpu==(i-1));     
    
    elem = intelem;
    for j=1:overlappinglevel
        % TO DO: C++ version of node2elem
        elem = node2elem(t(:,elem)+1,re,ce);
    end
    extelem = setdiff(elem,intelem);    
    extintelem = [intelem; extelem];            

    face = t2f(:,extintelem)+1;  
    face = unique(face(:));
    extface = setdiff(face,intface);        
    extintface = [intface; extface];

    intent = mkintent(facecon(:,extintface)+1,length(intface),face2cpu(extintface),i-1);         
    ent2cpu(intent) = i-1;                
    
    nent = nent + length(intent);    
end

if nent~=entmax
    error('EDG partition is incorrect');
end
if min(ent2cpu(:))<0
    error('EDG partition is incorrect');
end



