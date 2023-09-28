function [ds, dst, da, dat] = createGeoMatsTE(msh)
% CREATE_GEO_MATS_2D Creates the geometrical matrices in the TE case.
%
% Inputs:
%   msh - Mesh struct.
%
%  Outputs:
%   ds  - Primary edge matrix of size (np,np).
%   dst - Dual edge matrix of size (2np,2np).
%   da  - Primary area matrix of size (2np,2np).
%   dat - Dual area matrix of size (np,np).

    nx = msh.nx;
    ny = msh.ny;
    nz = 1;
    np = msh.np;
    xm = msh.xmesh;
    ym = msh.ymesh;

    % Primäres Gitter
    dx = [xm(2:end)-xm(1:(end-1)), 0];
    dy = [ym(2:end)-ym(1:(end-1)), 0];
    dz = 1;

    dsx = repmat(dx, 1, ny*nz);
    dsy = reshape(repmat(dy, nx, nz), 1, np);
    dsz = reshape(repmat(dz, nx*ny, 1), 1, np);

    ds = spdiags(dsz', 0, np, np);
    da = spdiags([dsy.*dsz, dsz.*dsx]', 0, 2*np, 2*np);

    % Duales Gitter
    dxt = [dx(1)/2, (dx(1:(end-2))+dx(2:(end-1)))/2, dx(end-1)/2];
    dyt = [dy(1)/2, (dy(1:(end-2))+dy(2:(end-1)))/2, dy(end-1)/2];

    dstx = repmat(dxt, 1, ny*nz);
    dsty = reshape(repmat(dyt, nx, nz), 1, np);

    dst = spdiags([dstx, dsty]', 0, 2*np, 2*np);
    dat = spdiags((dstx.*dsty)', 0, np, np);

end
