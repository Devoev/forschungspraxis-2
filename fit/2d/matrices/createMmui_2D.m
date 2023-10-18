function mmui = createMmui_2D(msh, ds, dst, da, rel_mui, mu0i)
% createMmui_2D creates the inverse permeability matrix for 2D mesh
%
% Input
% msh               -2D mesh object
% ds                -primary edge matrix
% dst               -dual edge matrix
% da                -primary face matrix
% rel_mui           -vector with distribution of relative inverse permeability
% mu0i              -inverse permeability of vacuum
%
% Output
% mmui              -inverse permeability matrix


% Get basic mesh parameters
np = msh.np;
nx = msh.nx;
ny = msh.ny;

% Include permeability of vacuum
mui = rel_mui * mu0i;

% Averaging permeability according to distribution of relative permeability
zeromat = sparse(np,np);
diagos = ones(np,2);
avex = spdiags(diagos, [0 -1], np, np);
avey = spdiags(diagos, [0 -nx], np, np);
avez = spdiags(diagos, [0 -nx*ny], np, np);

ave = [avex zeromat zeromat;
        zeromat avey zeromat;
        zeromat zeromat avez];

mui = reshape(mui,np,1);
mui = [ mui;mui;mui ];

% Calculate full vector for inverse permeability
dmui =  spdiags(0.5 * nullInv(dst) * ave * ds * mui, 0, 3*np, 3*np) .* [ones(2*np,1); 2*ones(np,1)];

% Calculate permeability matrix
mmui = dst * dmui * nullInv(da);
    
end
