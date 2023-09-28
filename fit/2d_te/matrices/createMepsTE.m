function meps = createMepsTE(msh, ds, da, dat, eps, bc)
% CREATE_MEPS_TE Creates the permittivity material matrix in the 2D TE case.
%
% Inputs:
%   msh  - Mesh struct.
%   ds   - Primary edge matrix.
%   da   - Primary area matrix.
%   dat  - Dual area matrix.
%   eps  - Permittivity.
%   bc   - Boundary conditions. TODO
%
% Outputs:
%   meps - Permittivity material matrix of size (np,np).

    if nargin < 6
        bc = [ 0 0 0 0 ];
    end
    if nargin < 5
        warning('Missing input parameters!')
    end

    nx = msh.nx;
    ny = msh.ny;
    np = msh.np;

    if numel(eps)==2

        % Homogener anisotroper Fall
        deps = spdiags([eps(1)*ones(np,1); eps(2)*ones(np,1)],0,np,np);

    elseif numel(eps)==1

        % Homogener isotroper Fall
        deps = eps*speye(np);

    elseif numel(eps)==np

        % Inhomogener Fall: Mittelung Permittivität
        zeromat = sparse(np,np);
        diagos = ones(np,4);
        avex = spdiags(diagos, [0 -nx -nx*ny -(nx+nx*ny)], np, np);
        avey = spdiags(diagos, [0 -nx*ny -1 -(nx*ny+1)], np, np);

        ave = [avex zeromat; zeromat avey];


	    eps = reshape(eps,np,1);
        eps = [eps;eps];

        deps = spdiags(0.25*nullInv(dat)*ave*da*eps,0,np,np);

    else
        warning('Nicht implementiert; falscher Übergabeparameter epsilon (Skalar oder Vektor mit Größe 2, np oder 2*np).')
    end

    % BC
    Mx = 1; My = nx;

    % BC in x direction
    if bc(1) || bc(2)
        indy = 1:ny;

        if bc(1)
            n=1+(1-1)*Mx+(indy-1)*My;
            deps(n,n) = 0;
        end
        if bc(2)
            n=1+(nx-1)*Mx+(indy-1)*My;
            deps(n,n) = 0;
        end
    end

    % BC in y direction
    if bc(3) || bc(4)
        indx = 1:nx;

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