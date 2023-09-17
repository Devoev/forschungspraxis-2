% Erzeugt die geometrischen Matrizen für ein kanonisches 
% kartesisches Gitter.
%
%   Eingabe
%   msh         Kanonisches kartesisches Gitter
%
%   Rückgabe
%   ds          Primäre Kanten-Matrix
%   dst         Duale Kanten-Matrix
%   da          Primäre Flächen-Matrix
%   dat         Duale Flächen-Matrix

function [ds, dst, da, dat] = createGeoMats(msh)

    nx = msh.nx;
    ny = msh.ny;
    nz = msh.nz;
    np = msh.np;
    xm = msh.xmesh;
    ym = msh.ymesh;
    zm = msh.zmesh;
        
    % Primäres Gitter
    dx = [xm(2:end)-xm(1:(end-1)), 0];
    dy = [ym(2:end)-ym(1:(end-1)), 0];
    dz = [zm(2:end)-zm(1:(end-1)), 0];
    
    dsx = repmat(dx, 1, ny*nz);
    dsy = reshape(repmat(dy, nx, nz), 1, np);
    dsz = reshape(repmat(dz, nx*ny, 1), 1, np);
    
    ds = spdiags([dsx, dsy, dsz]', 0, 3*np, 3*np);
    da = spdiags([dsy.*dsz, dsz.*dsx, dsx.*dsy]', 0, 3*np, 3*np);
        
    % Duales Gitter
    dxt = [dx(1)/2, (dx(1:(end-2))+dx(2:(end-1)))/2, dx(end-1)/2];
    dyt = [dy(1)/2, (dy(1:(end-2))+dy(2:(end-1)))/2, dy(end-1)/2];
    dzt = [dz(1)/2, (dz(1:(end-2))+dz(2:(end-1)))/2, dz(end-1)/2];
    
    dstx = repmat(dxt, 1, ny*nz);
    dsty = reshape(repmat(dyt, nx, nz), 1, np);
    dstz = reshape(repmat(dzt, nx*ny, 1), 1, np);
    
    dst = spdiags([dstx, dsty, dstz]', 0, 3*np, 3*np);
    dat = spdiags([dsty.*dstz, dstz.*dstx, dstx.*dsty]', 0, 3*np, 3*np);

end
