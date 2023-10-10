function [epart, npart] = elconpartition(elcon,np)

elcon = elcon';

numTotalEntities = max(elcon(:));
numRealEntities = length(unique(elcon(:)));
if numRealEntities ~= numTotalEntities
    error('There are ghost (missing) entities in the mesh.');
end

% number of elements
nt = size(elcon,1);

% current directory
cdir = pwd;

% move to directory that contains metis programs
if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end
ps=strcat(pwd,sslash);
is=find(ps==sslash);
up=0;
while ~(strcmp(ps(is(end-up-1)+1:is(end-up)-1),'hdgv1.0') || strcmp(ps(is(end-up-1)+1:is(end-up)-1),'HDGv1.0'))
    up = up+1;
end
cd(strcat(ps(1:is(end-up)),'metis'));

% generate a temporary file to be used in metis
dlmwrite('temp.txt', nt, 'delimiter', ' ','precision',10);
dlmwrite('temp.txt', elcon, '-append', 'delimiter', ' ','precision',10);

% call mpmetis
str = ['!./mpmetis temp.txt ' num2str(np)];
eval(str);

% get mesh partitioning data
str = ['temp.txt.epart.' num2str(np)];
epart = textread(str,'%d');

% get node partitioning data
str = ['temp.txt.npart.' num2str(np)];
npart = textread(str,'%d');

% remove files
%delete('temp.txt');
str = ['temp.txt.epart.' num2str(np)];
delete(str);
str = ['temp.txt.npart.' num2str(np)];
delete(str);

% move back to current directory
cd(cdir);


