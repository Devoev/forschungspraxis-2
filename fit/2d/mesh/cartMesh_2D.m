function [ msh ] = cartMesh_2D(xmesh, ymesh)
% cartMesh_2D generates the basic mesh object for a 2D domain
%
% Input
% xmesh             -1D vector with x-coordinates of the mesh
% ymesh             -1D vector with y-coordinates of the mesh
%
% Output
% msh               -Mesh object (struct) including:
%                    nx = number of points in x-direction
%                    ny = number of points in y-direction
%                    np = total number of points
%                    xmesh = same as input
%                    ymesh = same as input
%                    Mx = increment in x-direction of canonical indices
%                    My = increment in y-direction of canonical indices
%                    lz = length  of all edges in z-direction


% Determine nx, ny and np as well as Mx and My
nx = length(xmesh);
ny = length(ymesh);
np = nx*ny;
Mx = 1;
My = nx;

% Determine length lz of edges in z-direction
lz = min(min(diff(xmesh)), min(diff(ymesh)));

% Assign elements to struct msh
msh.nx = nx;
msh.ny = ny;
msh.np = np;
msh.Mx = Mx;
msh.My = My;
msh.xmesh = xmesh;
msh.ymesh = ymesh;
msh.lz = lz;

end
