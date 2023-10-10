function [t2f,t2t,f,ne,nf,nfe,nve,nvf] = mkt2ffast(t,elementtype,dim, check)
%MKT2T Compute Element to Element Connectivity.
%   T2T=MKT2T(T,MESHTYPE)
%
%      T:         Triangle indices (NT,3)
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%
%      T2T:       Triangle to Trangle Connectivity (NT,3)
%                 T2T(IT,IN) is the trangle that shares an edge
%                 with triangle IT and does nont contain node T(IT,IN).
%                 When an element is adjacent to the boundary the
%                 corresponding entry in T2T is set to zero
%

if nargin<=3 
    check = 0;
end

elem = unique(elementtype);
nelem = length(elem);
nfe = zeros(nelem);
nve = nfe;
nvf = nfe;
for i=1:nelem
    elemtype = elem(i);    
    switch dim
        case 1
            nfe(i) = 2;
            nve(i) = 2;
            nvf(i) = 1;
        case 2
            if elemtype==0
                nfe(i) = 3;
                nve(i) = 3;
                nvf(i) = 2;
            elseif elemtype==1
                nfe(i) = 4;
                nve(i) = 4;
                nvf(i) = 2;
            end
        case 3            
            if elemtype==0
                nfe(i) = 4;
                nve(i) = 4;
                nvf(i) = 3;
            elseif elemtype==1
                nfe(i) = 6;
                nve(i) = 8;
                nvf(i) = 4;
            elseif elemtype==2
                nfe(i) = 5;
                nve(i) = 6;
                nvf(i) = 4;
            elseif elemtype==3
                nfe(i) = 5;
                nve(i) = 5;
                nvf(i) = 4;
            end
        otherwise
            error('Only can handle dim=1, dim=2 or dim=3');
    end    
end

ne = numel(t)/max(nve);
sz = size(t);
if sz(1)==ne
    [t2t,t2f,f] = mkt2fc(t', elementtype, [dim ne max(nve) max(nfe) max(nvf)]);
else
    [t2t,t2f,f] = mkt2fc(t, elementtype, [dim ne max(nve) max(nfe) max(nvf)]);
end
f(:,f(1,:)==-1) = [];
nf = size(f,2);

if check == 1
    if sz(1)==ne
        [fm,t2fm,t2tm] = mkt2f(t+1,elemtype);
    else
        [fm,t2fm,t2tm] = mkt2f(t'+1,elemtype);
    end    
    if max(max(abs(t2t+1-double(t2tm)')))>1e-12
        error('t2t is incorrect.');
    end
    if max(max(abs(t2f+1-t2fm')))>1e-12
        error('t2f is incorrect.');
    end
    if max(max(abs(f(2:end,:)+1-fm')))>1e-12
        error('f is incorrect.');
    end
end


function [f,t2f,t2t] = mkt2f(t,elemtype)

if nargin<2, elemtype=0; end

[nt,npv]=size(t);
if npv==2 % 1D
    dim=1;
    nfv=2;
    npf=1;
else
    if elemtype==0 % tri/tet elements
        dim=size(t,2)-1;        
        nfv=dim+1;
        npf = dim;
    else % quad/hex elements
        dim=log2(size(t,2));        
        nfv=2*dim;
        npf=2*(dim-1);
    end
end

switch dim
    case 1
        face=[1;2];
    case 2
        if elemtype==0
            face=[[2,3];[3,1];[1,2]];
        elseif elemtype==1
            face=[[1,2];[2,3];[3,4];[4,1]];
        end
    case 3
        if elemtype==0
            face=[[2,3,4];[1,4,3];[1,2,4];[1,3,2]];
        elseif elemtype==1
            face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]];
        end
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end


t2t = mkt2t(t,elemtype);
nb = sum(sum(t2t <= 0));
f = zeros((nfv*nt+nb)/2,npf+2);
t2f = zeros(nt,nfv);
jf = 0;
for i=1:nt
    for j=1:nfv
        if t2t(i,j) > i || t2t(i,j) <=0
            ie = t2t(i,j);
            jf = jf + 1;
            
            f(jf,1:npf) = t(i,face(j,:));
            f(jf,npf+1) = i;
            f(jf,npf+2) = ie;
            t2f(i,j) = jf;
            
            if ie > 0
                %k = sum(reshape(t(ie,face),[nfv npf]),2)-sum(f(jf,1:npf))==0;                                                            
                %t2f(ie,k) = jf;                
                a = sort(reshape(t(ie,face),[nfv npf]),2);
                b = sort(f(jf,1:npf));                
                k = sum(abs(a-repmat(b,[nfv 1])),2)==0;
                t2f(ie,k) = jf;                                
            end                        
        end
    end
end


