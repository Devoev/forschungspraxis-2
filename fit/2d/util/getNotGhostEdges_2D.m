function idxNotGhostEdges = getNotGhostEdges_2D(msh)
% getNotGhostEdges_2D gets the indices of all ghost edges on the calculation
% domain and saves them in the vector idxGhostEdges
%
% Inputs:
%   msh                 - Mesh struct
%
% Outputs:
%   idxNotGhostEdges    - Vector, containing the indices of all 
%                         edges relevant for the calculation


% Get basic mesh parameters
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
n_ymax = 1+(indx-1)*Mx+(ny-1)*My;

% return all indices of ghost edges in one column vector
idxGhostEdges = [n_xmax,np+n_ymax]';

% Determine indices of calc relevant edges
idxNotGhostEdges = setdiff((1:3*np), idxGhostEdges);
idxNotGhostEdges = idxNotGhostEdges';


end
