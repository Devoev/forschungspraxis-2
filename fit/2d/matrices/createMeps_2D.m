
function meps = createMeps_2D(msh, ds, da, dat, rel_eps, eps0)
%% Description
%
% Create permittivity matrix for 2D mesh
%
% Input
% msh               2D mesh object
% ds                primary edge matrix
% da                primary face matrix
% dat               dual face matrix
% rel_eps           vector with distribution of relative permeability
% eps0              permeability of vacuum
%
% Output
% meps              permittivity matrix


%% Function definition

    np = msh.np;
    nx = msh.nx;
    ny = msh.ny;

    % Combine relative and vacuum permeability
    epsilon = rel_eps * eps0;

    % Calculate vector for avaraged epsilon
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

    deps = spdiags(0.25*nullInv(dat)*ave*da*epsilon,0,3*np,3*np) .* [2*ones(2*np,1); ones(np,1)];

    % Calculate permittivity matrix
    meps = dat * deps * nullInv(ds);
    
end
