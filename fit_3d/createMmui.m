%Erzeugt die Materialmatrix der inversen Permeabilität.
%
%   Eingabe
%   msh         Kanonisches kartesisches Gitter
%   ds          Primäre Kanten-Matrix
%   dst         Duale Kanten-Matrix
%   da          Primäre Flächen-Matrix
%   mui         Inverse Permeabilität
%   bc          bc=0 für magnetische Randbedingungen
%               bc=1 für elektrische Randbedingungen
%
%   Rückgabe
%   mmui        Materialmatrix der inversen Permeabilitäten

function mmui = createMmui(msh, ds, dst, da, mui, bc)
    
    if nargin < 6
        bc = [ 0 0 0 0 0 0 ];
    end
    if nargin < 5
        warning('Zu wenig Eingabeparameter.')
    end
    
    nx = msh.nx;    ny = msh.ny;    nz = msh.nz;
    np = msh.np; 
   
    if numel(mui)==3
        
        % Homogener anisotroper Fall
        dmui = [mui(1)*ones(np,1); mui(2)*ones(np,1); mui(3)*ones(np,1)];
        
    elseif numel(mui)==1
        
        % Homogener isotroper Fall
        dmui = mui*ones(3*np,1);
        
    else

        % Inhomogener Fall: Mittelung Permeabilitäten
        zeromat = sparse(np,np);
        diagos = ones(np,2);
        avex = spdiags(diagos, [0 -1], np, np);
        avey = spdiags(diagos, [0 -nx], np, np);
        avez = spdiags(diagos, [0 -nx*ny], np, np);

        ave = [avex zeromat zeromat;
                zeromat avey zeromat;
                zeromat zeromat avez];
        
        if numel(mui) == np
            mui = reshape(mui,np,1);
             mui = [ mui;mui;mui ];
        elseif numel(mui) ~= 3*np
            warning('Nicht implementiert; falscher Übergabeparameter mui (Skalar oder Vektor mit Größe 3, np oder 3*np).')
        end

        dmui =  0.5 * nullInv(dst) * ave * ds * mui;
        
    end
    
    % Randbedingungen einarbeiten
    Mx = 1; My = nx; Mz = nx*ny;
    
    if bc(1) || bc(2)
        indy = repmat(1:ny,1,nz);
        indz = reshape(repmat(1:nz,ny,1),1,ny*nz);
        if bc(1)
            n=1+(1-1)*Mx+(indy-1)*My+(indz-1)*Mz;
            dmui(n) = 0;
        end
        if bc(2)
            n=1+(nx-1)*Mx+(indy-1)*My+(indz-1)*Mz;
            dmui(n) = 0;
        end
    end
    if bc(3) || bc(4)
        indx = repmat(1:nx,1,nz);
        indz = reshape(repmat(1:nz,nx,1),1,nx*nz);
        if bc(3)
            n = 1+(indx-1)*Mx+(1-1)*My+(indz-1)*Mz; 
            dmui(n+np) = 0;
        end
        if bc(4)
            n = 1+(indx-1)*Mx+(ny-1)*My+(indz-1)*Mz;
            dmui(n+np) = 0;
        end
    end
    if bc(5) || bc(6)
        indx = repmat(1:nx,1,ny);
        indy = reshape(repmat(1:ny,nx,1),1,nx*ny);
        if bc(5)
            n=1+(indx-1)*Mx+(indy-1)*My+(1-1)*Mz;
            dmui(n+2*np) = 0;
        end
        if bc(6)
            n=1+(indx-1)*Mx+(indy-1)*My+(nz-1)*Mz;
            dmui(n+2*np) = 0;
        end
    end
    
    % geometrische Matritzen ranmuliplizieren
    mmui = dst * spdiags(dmui,0,3*np,3*np) * nullInv(da);
    
end
