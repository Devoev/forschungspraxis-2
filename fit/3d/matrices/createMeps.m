% Erzeugt die geometrischen Matrizen für ein kanonisches 
% kartesisches Gitter.
%
%   Eingabe
%   msh         Kanonisches kartesisches Gitter
%   ds          Primäre Kanten-Matrix
%   da          Primäre Flächen-Matrix
%   dat         Duale Flächen-Matrix
%   epsilon     Relative Permittivität als Skalar oder Vektor an. Vektor kann
%               folgende Größe haben
%               1 x 3 haben für homogen, aber isotrop
%               1 x Np für inhomogen, aber anisotrop
%               1 x 3*Np für inhomogen und anistrop.
%   bc          Randbedingungen (0: magnetisch, 1: elektrisch)
%
%   Rückgabe
%   meps        Materialmatrix der Permittivitäten

function meps = createMeps(msh, ds, da, dat, epsilon, bc)
    
    if nargin < 6
        bc = [ 0 0 0 0 0 0 ];
    end
    if nargin < 5
        warning('Zu wenig Eingabeparameter.')
    end
    
    nx = msh.nx;
    ny = msh.ny;
    nz = msh.nz;
    np = msh.np;
   
    if numel(epsilon)==3
        
        % Homogener anisotroper Fall
        deps = spdiags([epsilon(1)*ones(np,1); epsilon(2)*ones(np,1); epsilon(3)*ones(np,1)],0,3*np,3*np);

    elseif numel(epsilon)==1
        
        % Homogener isotroper Fall
        deps = epsilon*speye(3*np);
        
    elseif numel(epsilon)==np
    
        % Inhomogener Fall: Mittelung Permittivität
        zeromat = sparse(np,np);
        diagos = ones(np,4);
        avex = spdiags(diagos, [0 -nx -nx*ny -(nx+nx*ny)], np, np);
        avey = spdiags(diagos, [0 -nx*ny -1 -(nx*ny+1)], np, np);
        avez = spdiags(diagos, [0 -1 -nx -(1+nx)], np, np);

        ave = [avex zeromat zeromat;
                zeromat avey zeromat;
                zeromat zeromat avez];
        
	epsilon = reshape(epsilon,np,1);
        epsilon = [ epsilon;epsilon;epsilon ];

        deps = spdiags(0.25*nullInv(dat)*ave*da*epsilon,0,3*np,3*np);

    else
        warning('Nicht implementiert; falscher Übergabeparameter epsilon (Skalar oder Vektor mit Größe 3, np oder 3*np).')
    end
    
    % Randbedingungen einarbeiten
    Mx = 1; My = nx; Mz = nx*ny;
    
    if bc(1) || bc(2)
        indy = repmat(1:ny,1,nz);
        indz = reshape(repmat(1:nz,ny,1),1,ny*nz);
        if bc(1)
            n=1+(1-1)*Mx+(indy-1)*My+(indz-1)*Mz;
            deps(n+np,n+np) = 0;
            deps(n+2*np,n+2*np) = 0;
        end
        if bc(2)
            n=1+(nx-1)*Mx+(indy-1)*My+(indz-1)*Mz;
            deps(n+np,n+np) = 0;
            deps(n+2*np,n+2*np) = 0;
        end
    end
    if bc(3) || bc(4)
        indx = repmat(1:nx,1,nz);
        indz = reshape(repmat(1:nz,nx,1),1,nx*nz);
        if bc(3)
            n = 1+(indx-1)*Mx+(1-1)*My+(indz-1)*Mz; 
            deps(n,n) = 0;
            deps(n+2*np,n+2*np) = 0;
        end
        if bc(4)
            n = 1+(indx-1)*Mx+(ny-1)*My+(indz-1)*Mz;
            deps(n,n) = 0;
            deps(n+2*np,n+2*np) = 0;
        end
    end
    if bc(5) || bc(6)
        indx = repmat(1:nx,1,ny);
        indy = reshape(repmat(1:ny,nx,1),1,nx*ny);
        if bc(5)
            n=1+(indx-1)*Mx+(indy-1)*My+(1-1)*Mz;
            deps(n,n) = 0;
            deps(n+np,n+np) = 0;
        end
        if bc(6)
            n=1+(indx-1)*Mx+(indy-1)*My+(nz-1)*Mz;
            deps(n,n) = 0;
            deps(n+np,n+np) = 0;
        end
    end
    
    
    meps = dat*deps*nullInv(ds);
    
end