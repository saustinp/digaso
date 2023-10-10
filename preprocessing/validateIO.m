filename = [app.filein 'app'];
fileID = fopen([filename,'.bin'],'r');
a1 = fread(fileID,'double');
fclose(fileID);

filename = [app.fileout 'app'];
fileID = fopen([filename,'.bin'],'r');
a2 = fread(fileID,'double');
fclose(fileID);

e = max(abs(a1(:)-a2(:)));
if e>1e-12
    error('app struct is incorrect.');
else
    disp('app struct is correct.');
end

filename = [app.filein 'master'];
fileID = fopen([filename,'.bin'],'r');
a1 = fread(fileID,'double');
fclose(fileID);

filename = [app.fileout 'master'];
fileID = fopen([filename,'.bin'],'r');
a2 = fread(fileID,'double');
fclose(fileID);

e = max(abs(a1(:)-a2(:)));
if e>1e-12
    error('master struct is incorrect.');
else
    disp('master struct is correct.');    
end

for i = 1:app.nfile
    filename = [app.filein 'mesh' num2str(i)];
    fileID = fopen([filename,'.bin'],'r');
    a1 = fread(fileID,'double');
    fclose(fileID);

    filename = [app.fileout 'mesh' num2str(i)];
    fileID = fopen([filename,'.bin'],'r');
    a2 = fread(fileID,'double');
    fclose(fileID);
        
    e = max(abs(a1(:)-a2(:)));
    if e>1e-12
        error('mesh struct is incorrect.');        
    else
        disp('mesh struct is correct.');    
    end        
end

for i = 1:app.nfile
    filename = [app.filein 'sol' num2str(i)];
    fileID = fopen([filename,'.bin'],'r');
    a1 = fread(fileID,'double');
    fclose(fileID);

    filename = [app.fileout 'sol' num2str(i)];
    fileID = fopen([filename,'.bin'],'r');
    a2 = fread(fileID,'double');
    fclose(fileID);
        
    e = max(abs(a1(:)-a2(:)));
    if e>1e-12
        error('sol struct is incorrect.');        
    else
        disp('sol struct is correct.');    
    end        
end


