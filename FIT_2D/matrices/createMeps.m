%% Description
%
% Create permittivity matrix for 2D mesh
%
% Input
% msh               2D mesh object
% ds                primary edge matrix
% dat               dual face matrix
% epsilon           value for relative epsilon on domain
%
% Output
% meps              permittivity matrix


%% Function definition
function meps = createMeps(msh, ds, dat, epsilon)

    % Include permittivity of vacuum
    epsilon = epsilon * 8.854187e-12;

    % Homogenious isotrope case
    deps = epsilon * speye(3 * msh.np);
   
    % Calculate permittivity matrix
    meps = dat * deps * nullInv(ds);
    
end
