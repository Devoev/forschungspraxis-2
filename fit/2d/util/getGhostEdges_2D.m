function idxGhostEdges = getGhostEdges_2D(msh)
% getGhostEdges_2D gets the indices of all ghost edges on the calculation
% domain and saves them in the vector idxGhostEdges
%
% Inputs:
%   msh             - Mesh struct
%
% Outputs:
%   idxGhostEdges   - Vector, containing the indices of all ghost edges


% Get basic mesh parameters
np = msh.np;
nx = msh.nx;
ny = msh.ny;
Mx = msh.Mx;
My = msh.My;

% Calculate the indices of ghost edges in x-direction
indy = repmat(1:ny,1,1);
n_xmax = 1+(nx-1)*Mx+(indy-1)*My;

% Calculate the indices of ghost edges in y-direction
indx = repmat(1:nx,1,1);
n_ymax = 1+(indx-1)*Mx+(ny-1)*My;

% Return all indices of ghost edges in one column vector
idxGhostEdges = [n_xmax,np+n_ymax]';


end
