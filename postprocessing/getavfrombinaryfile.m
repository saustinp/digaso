
function av = getavfrombinaryfile(foldername, filename, avfilename, nproc, timeStepNo, npv, ncav)

if nproc==1
    fileID = fopen([foldername '/' avfilename '_t' num2str(timeStepNo) '.bin'],'r');
    av = reshape(fread(fileID,ndims(1),'double'),npv,[]);
    fclose(fileID);
else
    elempart = [];
    nelem = zeros(nproc,1);
    for i = 1:nproc
        fileID = fopen([foldername '/' filename 'elempart_np' num2str(i-1) '.bin'],'r');
        elemparti = fread(fileID,'int');
        elemparti = elemparti(1:2:end-1);
        nelem(i) = length(elemparti);
        elempart = [elempart; elemparti];
        fclose(fileID);
    end
    elempart = elempart + 1;
    
    n1 = cumsum([0; nelem]); 
    
    av = zeros(npv,ncav,n1(end));      
    for i = 1:nproc        
        fullFileName = [foldername  '/' avfilename '_t' num2str(timeStepNo) '_np' num2str(i-1) '.bin'];
        fileID = fopen(fullFileName,'r');    
        
        ne = nelem(i);                        
        av(:,:,elempart((n1(i)+1):(n1(i)+ne))) = reshape(fread(fileID,'double'),[npv ncav ne]);
        
        fclose(fileID);
    end
end
