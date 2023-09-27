
function idxNotGhostEdges = getNotGhostEdges_2D(msh)
%% Description
%
% Function to get indices of calc relevant edges (not ghost edges)
%
% Input
% msh               2D mesh object
%
% Output
% idxNotGhostEdges  Indices of calc relevant edges    


%% Function definition

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
    idxNotGhostEdges = linspace(1, 3*np, 3*np);
    idxNotGhostEdges(idxGhostEdges) = 0;
    idxNotGhostEdges = idxNotGhostEdges';

end
