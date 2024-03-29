function mmui = createMmuiTE(msh, ds, dst, da, mui, bc)
% CREATE_MMUI_TE Creates the reluctivity material matrix in the 2D TE case.
%
% Inputs:
%   msh  - Mesh struct.
%   ds   - Primary edge matrix.
%   dst  - Dual edge matrix.
%   da   - Primary area matrix.
%   mui  - Inverse permeability.
%   bc   - Boundary conditions. TODO
%
% Outputs:
%   mmui - Reluctivity material matrix of size (2np,2np).
    
    if nargin < 6
        bc = [ 0 0 0 0 ];
    end
    if nargin < 5
        warning('Missing input parameters!')
    end
    
    nx = msh.nx;
    ny = msh.ny;
    np = msh.np; 
   
    if numel(mui)==2
        
        % Homogener anisotroper Fall
        dmui = [mui(1)*ones(np,1); mui(2)*ones(np,1)];
        
    elseif numel(mui)==1
        
        % Homogener isotroper Fall
        dmui = mui*ones(2*np,1);
        
    else

        % Inhomogener Fall: Mittelung Permeabilitäten
        zeromat = sparse(np,np);
        diagos = ones(np,2);
        avex = spdiags(diagos, [0 -1], np, np);
        avey = spdiags(diagos, [0 -nx], np, np);
    

        ave = [avex zeromat; zeromat avey ];

            
        
        if numel(mui) == np
            mui = reshape(mui,np,1);
             mui = [ mui;mui];
        elseif numel(mui) ~= 2*np
            warning('Nicht implementiert; falscher Übergabeparameter mui (Skalar oder Vektor mit Größe 2, np oder 2*np).')
        end

        dmui =  0.5 * nullInv(dst) * ave * ds * mui;
        
    end
    
    % Randbedingungen einarbeiten
    Mx = 1; My = nx; 

    % BC in x direction
    if bc(1) || bc(2)
        indy = 1:ny;

        if bc(1)
            n=1+(1-1)*Mx+(indy-1)*My;
            dmui(n) = 0;
        end
        if bc(2)
            n=1+(nx-1)*Mx+(indy-1)*My;
            dmui(n) = 0;
        end
    end

    % BC in y direction
    if bc(3) || bc(4)
        indx = 1:nx;

        if bc(3)
            n = 1+(indx-1)*Mx+(1-1)*My; 
            dmui(n+np) = 0;
        end
        if bc(4)
            n = 1+(indx-1)*Mx+(ny-1)*My;
            dmui(n+np) = 0;
        end
    end
    
    
    % geometrische Matritzen ranmuliplizieren
    mmui = dst * spdiags(dmui,0,2*np,2*np) * nullInv(da);

end
