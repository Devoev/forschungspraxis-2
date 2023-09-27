%% Description
%
% Create geometry matrices for 2D mesh
%
% Input
% msh               2D mesh object
%
% Output
% ds                matrix for primary edges
% dst               matrix for dual edges
% da                matrix for primary faces
% dat               matrix for dual faces


%% Function definition
function [ds, dst, da, dat] = createGeoMats(msh)

    % Extraction of input parameters
    nx = msh.nx;
    ny = msh.ny;
    np = msh.np;
    xm = msh.xmesh;
    ym = msh.ymesh;
        
    % Primary mesh
    dx = [xm(2:end)-xm(1:(end-1)), 0];
    dy = [ym(2:end)-ym(1:(end-1)), 0];
    
    dsx = repmat(dx, 1, ny);
    dsy = reshape(repmat(dy, nx, 1), 1, np);
        
    ds = spdiags([dsx, dsy, ones(1, np)]', 0, 3*np, 3*np);
    da = spdiags([dsy, dsx, dsx.*dsy]', 0, 3*np, 3*np);
        
    % Dual mesh
    dxt = [dx(1)/2, (dx(1:(end-2))+dx(2:(end-1)))/2, dx(end-1)/2];
    dyt = [dy(1)/2, (dy(1:(end-2))+dy(2:(end-1)))/2, dy(end-1)/2];
    
    dstx = repmat(dxt, 1, ny);
    dsty = reshape(repmat(dyt, nx, 1), 1, np);
    
    dst = spdiags([dstx, dsty, ones(1, np)]', 0, 3*np, 3*np);
    dat = spdiags([dsty, dstx, dstx.*dsty]', 0, 3*np, 3*np);

end
