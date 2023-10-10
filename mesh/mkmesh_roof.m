function mesh = mkmesh_roof(m,n,o,porder,a,b,c,elemtype,nodetype)
%MKMESH_CUBE Creates 3D mesh data structure for unit square.
%   MESH=MKMESH_CUBE(M,N,PORDER,PARITY)
%
%      MESH:      Mesh structure
%      M:         Number of points in the x direction 
%      N:         Number of points in the y direction
%      O:         Number of points in the z direction
%      PORDER:    Polynomial Order of Approximation (default=1)
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%      NODETYPE:  Flag determining node distribution 
%                 Flag = 0 uniform distribution (default)
%                 Flag = 1 nonuniform distribution

%   See also: CUBEMESH, MKMESH
%

if nargin<1, m=2; end
if nargin<2, n=m; end
if nargin<3, o=n; end
if nargin<4, porder = 1; end
if nargin<5, a=1;      end
if nargin<6, b=1;      end
if nargin<7, c=1;      end
if nargin<8, elemtype=0; end
if nargin<9, nodetype=0; end


if m < 2 || n < 2 || o < 2,
    error('At least m=2, n=2, o=2 needed.');
end

[p,t] = cubemesh(m,n,o,elemtype);
p(:,1) = a*p(:,1);
p(:,2) = b*p(:,2);
p(:,3) = c*p(:,3);

bndexpr = {'all(p(:,1)<1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)<1e-3)','all(p(:,2)>max(p0(:,2))-1e-3)', ...
           'all(p(:,3)<1e-3)','all(p(:,3)>max(p0(:,3))-1e-3)'};     
       
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

alfa = 40*pi/180;
R0 = 25;
R = R0+mesh.p(:,3);
%x = R.*sin(alfa).*mesh.p(:,1);
x = (R0.*sin(alfa)+mesh.p(:,3)*sin(alfa)).*mesh.p(:,1)-mesh.p(:,3)*sin(alfa);
z = -R0.*cos(alfa)+sqrt(R.*R-(x-R0.*sin(alfa)).^2);
mesh.p(:,3) = z;
mesh.p(:,1) = x;

R = R0+mesh.dgnodes(:,3,:);
%x = R.*sin(alfa).*mesh.dgnodes(:,1,:);
x = (R0.*sin(alfa)+mesh.dgnodes(:,3,:)*sin(alfa)).*mesh.dgnodes(:,1,:)-mesh.dgnodes(:,3,:)*sin(alfa);
z = -R0.*cos(alfa)+sqrt(R.*R-(x-R0.*sin(alfa)).^2);
mesh.dgnodes(:,3,:) = z;
mesh.dgnodes(:,1,:) = x;

% alfa = pi/6;
% R = 1;
% x = mesh.p(:,1);
% z = -R*cos(alfa)+sqrt(R*R-(x-R*sin(alfa)).^2);
% mesh.p(:,3) = z+mesh.p(:,3);
% 
% x = mesh.dgnodes(:,1,:);
% z = -R*cos(alfa)+sqrt(R*R-(x-R*sin(alfa)).^2);
% mesh.dgnodes(:,3,:) = z+mesh.dgnodes(:,3,:);
