
function [ msh ] = cartMesh_2D(xmesh, ymesh)
%% Description
%
% Creates a basic 2D mesh object
%
% Input
% xmesh             1D vector with x-coordinates of the mesh
% ymesh             1D vector with y-coordinates of the mesh
%
% Output
% msh               Mesh object (struct) including:
%                   nx = number of points in x-direction
%                   ny = number of points in y-direction
%                   np = total number of points
%                   xmesh = same as input
%                   ymesh = same as input
%                   Mx = increment in x-direction of canonical indices
%                   My = increment in y-direction of canonical indices


%% Function definition

    % Determine nx, ny and np as well as Mx and My
    nx = length(xmesh);
    ny = length(ymesh);
    np = nx*ny;
    Mx = 1;
    My = nx;
    
    % Assign elements to struct msh
    msh.nx = nx;
    msh.ny = ny;
    msh.np = np;
    msh.Mx = Mx;
    msh.My = My;
    msh.xmesh = xmesh;
    msh.ymesh = ymesh;

end
