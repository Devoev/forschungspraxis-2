function idxGhostEdges = getGhostEdges2D(msh)
% getGhostEdges(msh) returns the canonical indices of all ghost edges in
% the considered mesh
%   
% Input:
%   msh             msh struct as created by cartMesh2D
%
% Output:
%   idxGhostEdges   column vector that contains all ghost edges' indices
%
% See also cartMesh

np = msh.np;
nx = msh.nx;
ny = msh.ny;
Mx = msh.Mx;
My = msh.My;

% calculates indices for boundary xmax
indy = repmat(1:ny,1,1);
n_xmax = 1+(nx-1)*Mx+(indy-1)*My;

% calculates indices for boundary ymax
indx = repmat(1:nx,1,1);
n_ymax = np + 1+(indx-1)*Mx+(ny-1)*My;

% return all indices in one column vector
idxGhostEdges = [n_xmax,n_ymax]';
