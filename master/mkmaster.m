function master=mkmaster(mesh,pgauss)
%MKMASTER  Create master element structure
%    MASTER=MKMASTER(MESH)
%
%      MESH:      Mesh data structure
%      PGAUSS:    Degree of the polynomila to be integrated exactly
%                 (default: PGAUSS = 2*MESH.PORDER)
%

if nargin < 2 
    pgauss = max(2*mesh.porder,1); 
end

master.nd     = mesh.nd;       % problem dimension
master.porder = mesh.porder;   % polynomial degree
master.plocvl = mesh.plocal;   % nodal points on the master volume element
master.plocfc = mesh.plocfc;   % nodal points on the master face element
master.perm   = mesh.perm;     % locations of the nodal points on the faces of the master elements

dim      = mesh.nd;
elemtype = mesh.elemtype;      % element type flag: simplex or tensor
porder   = mesh.porder;

% Gauss points and weights on the master volume element  
[master.gpvl,master.gwvl] = gaussquad(pgauss,dim,elemtype);
% master.gpnvl = master.plocvl;
% master.gwnvl = (1/size(master.plocvl,1)) * ones(size(master.plocvl,1),1);

% shape functions and derivatives on the master volume element  
master.shapvl = mkshape(porder,master.plocvl,master.gpvl,elemtype);
% master.shapnvl = mkshape(porder,master.plocvl,master.plocvl,elemtype);

if dim>1
    % Gauss points and weights on the master face element  
    [master.gpfc,master.gwfc] = gaussquad(pgauss,dim-1,elemtype);
%     master.gpnfc = master.plocfc;
%     master.gwnfc = (1/size(master.plocfc,1)) * ones(size(master.plocfc,1),1);

    % shape functions and derivatives on the master face element  
    master.shapfc = mkshape(porder,master.plocfc,master.gpfc,elemtype);
%     master.shapnfc = mkshape(porder,master.plocfc,master.plocfc,elemtype);
else
    master.plocfc = 0;
    master.gpfc   = 0;
    master.gwfc   = 1;
    master.shapfc = 1;
%     master.gpnfc   = 0;
%     master.gwnfc   = 1;
%     master.shapnfc = 1;
end

% mass and convection matrices on the master element
master.mass = squeeze(master.shapvl(:,:,1))*diag(master.gwvl)*squeeze(master.shapvl(:,:,1))';
for ii=1:dim
    master.conv(:,:,ii) = squeeze(master.shapvl(:,:,1))*diag(master.gwvl)*squeeze(master.shapvl(:,:,ii+1))';
end   

master.ngv = size(master.gpvl,1);   % number of gasss points per element
master.ngf = size(master.gpfc,1);   % number of gasss points per face
master.npv = size(master.plocvl,1); % number of nodal points per element
master.npf = size(master.plocfc,1); % number of nodal points per face
