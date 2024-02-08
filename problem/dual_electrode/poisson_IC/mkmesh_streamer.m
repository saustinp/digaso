function mesh = mkmesh_rect(m,n,porder,parity,xrect,elemtype,nodetype)
%MKMESH_SQUARE Creates 2D mesh data structure for unit square.
%   MESH=MKMESH_SQUARE(M,N,PORDER,PARITY)
%
%      MESH:      Mesh structure
%      M:         Number of points in the horizaontal direction 
%      N:         Number of points in the vertical direction
%      PORDER:    Polynomial Order of Approximation (default=1)
%      PARITY:    Flag determining the the triangular pattern
%                 Flag = 0 (diagonals SW - NE) (default)
%                 Flag = 1 (diagonals NW - SE)
%
%   See also: SQUAREMESH, MKT2F, SETBNDNBRS, UNIFORMLOCALPNTS, CREATENODES
%

% Adapted from mkmesh_rect

% m is number of horizontal elements, n is vertical
parity = 0;
xmax = 125;
ymax = 125;
elemtype = 0;
nodetype = 1;

[p,t] = squaremesh(m,n,parity,elemtype);
p(:,1) = p(:,1).^2 * xmax;      % Squaring the x coord shifts the points closer to the y axis, refining there
p(:,2) = ymax*p(:,2);

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     

mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
