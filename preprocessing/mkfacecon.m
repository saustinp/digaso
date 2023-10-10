function facecon = mkfacecon(elcon,f,t2f)

[ne,nfe] = size(t2f);
npf = numel(elcon)/(nfe*ne); 
nf = size(f,1);

elcon = reshape(elcon,[npf,nfe ne]);
facecon = zeros(npf,nf);
for i = 1:nf
    fe = f(i,end-1:end); % neighboring elements of face i        
    if1 = t2f(fe(1),:)==i; % location of face i on element fe(1)    
    facecon(:,i) = elcon(:,if1,fe(1));    
end
