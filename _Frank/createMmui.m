%% Description
%
% Create inverse permeability matrix for 2D mesh
%
% Input
% msh               2D mesh object
% dst               dual edge matrix
% da                primary face matrix
% mui               value for inverse relative permeability on domain
%
% Output
% mmui              inverse permeability matrix


%% Function definition
function mmui = createMmui(msh, dst, da, mui)

    % Include permeability of vacuum
    mui = mui / (pi*4e-7);

    % Homogenious isotrope case
    dmui = mui * speye(3 * msh.np);

    % Calculate permeability matrix
    mmui = dst * dmui * nullInv(da);
    
end
