% Erzeugt die geometrischen Matrizen für ein kanonisches 
% kartesisches Gitter.
%
%   Eingabe
%   msh         Kanonisches kartesisches Gitter
%   ds          Primäre Kanten-Matrix
%   da          Primäre Flächen-Matrix
%   dat         Duale Flächen-Matrix
%   eps     Relative Permittivität als Skalar oder Vektor an. 
%   bc          Randbedingungen (0: magnetisch, 1: elektrisch)
%
%   Rückgabe
%   meps        Materialmatrix der Permittivitäten

function meps = create2DMeps(msh, ds, da, dat, eps, bc)
    
    if nargin < 6
        bc = [ 0 0 0 0 0 0 ];
    end
    if nargin < 5
        warning('Zu wenig Eingabeparameter.')
    end
    
    nx = msh.nx;
    ny = msh.ny;
    np = msh.np;
   
    if numel(eps)==2
        
        % Homogener anisotroper Fall
        deps = spdiags([eps(1)*ones(np,1); eps(2)*ones(np,1)],0,2*np,2*np);

    elseif numel(eps)==1
        
        % Homogener isotroper Fall
        deps = eps*speye(2*np);
        
    elseif numel(eps)==np
    
        % Inhomogener Fall: Mittelung Permittivität
        zeromat = sparse(np,np);
        diagos = ones(np,4);
        avex = spdiags(diagos, [0 -nx -nx*ny -(nx+nx*ny)], np, np);
        avey = spdiags(diagos, [0 -nx*ny -1 -(nx*ny+1)], np, np);
        
        ave = [avex zeromat; zeromat avey];
             
        
	    eps = reshape(eps,np,1);
        eps = [eps;eps];

        deps = spdiags(0.25*nullInv(dat)*ave*da*eps,0,2*np,2*np);

    else
        warning('Nicht implementiert; falscher Übergabeparameter epsilon (Skalar oder Vektor mit Größe 2, np oder 2*np).')
    end
    
    % Randbedingungen einarbeiten
    Mx = 1; My = nx; 
    
    if bc(1) || bc(2)
        indy = repmat(1:ny,1,1);
     
        if bc(1)
            n=1+(1-1)*Mx+(indy-1)*My;
            deps(n+np,n+np) = 0;
            
        end
        if bc(2)
            n=1+(nx-1)*Mx+(indy-1)*My;
            deps(n+np,n+np) = 0;
            
        end
    end
    if bc(3) || bc(4)
        indx = repmat(1:nx,1,1);
      
        if bc(3)
            n = 1+(indx-1)*Mx+(1-1)*My; 
            deps(n,n) = 0;
         
        end
        if bc(4)
            n = 1+(indx-1)*Mx+(ny-1)*My;
            deps(n,n) = 0;
          
        end
    end
    
    
    
    meps = dat*deps*nullInv(ds);
    
end