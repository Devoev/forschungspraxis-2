function ebow_abs = calc_abs_field(msh, ebow)
% CALC_ABS_FIELD Calculates the absolute value of the given field ebow.
%
% Inputs:
%   msh         - Mesh struct.
%   ebow        - Phasor of integrated electrical field. Vector of size (3*np).
%
% Outputs:
%   ebow_abs    - Absolute value of electrical field. Vector of size (np).

    % Indices of edges
    idx_x = 1:msh.np;
    idx_y = 1+msh.np:2*msh.np;
    idx_z = 1+2*msh.np:3*msh.np;

    % Calc absolute value
    ebow_abs = sqrt(abs(ebow(idx_x)).^2 + abs(ebow(idx_y)).^2 + abs(ebow(idx_z)).^2);   % TODO: Scale with da/ds

end