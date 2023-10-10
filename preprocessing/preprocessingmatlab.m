function [elcon,t2f] = preprocessingmatlab(p, t, elementtype, nproc)

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

if nproc == 1
    [t2f,t2t,f,ne,nf,nfe,nve,nvf] = mkt2fmatlab(t'+1,elementtype,dim);    
    t2f = t2f'-1; t2t = t2t'-1; f(:,2:end) = f(:,2:end)-1; f = f';         
    [elcon,facecon,edg] = mkelconcpp(p,t,f,t2f,elementtype,isEDGface,app.porder,dim,check);            
    f = setbndnbrs(p,f,bndexpr);    
    bf = reshape(f(end,t2f+1),size(t2f));       
    if isempty(periodicexpr)==0
        f1 = f(1,:); f = f(2:end,:);        
        [elcon,t2f,bf,f,t,isEDGface] = periodicmatlab(p,elcon+1,t2f'+1,bf+1,f'+1,t'+1,isEDGface,...
                                       app.porder(1),elementtype(1),master.perm,periodicexpr,app.hybrid);        
        elcon = elcon-1; t2f = t2f'-1; bf = bf-1; f = f'-1; t = t'-1; 
        f = [f1; f];
    end    
else
    
end

