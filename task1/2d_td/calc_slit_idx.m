function idx = calc_slit_idx(msh, d, delta, use_y_symmetry, polarisation)
% CALC_SLIT_IDX Calculates the indices for the the 2 slits of the double slit.
%
% Inputs:
%   msh             - Mesh.
%   d               - Slit distance.
%   delta           - Slit width.
%   use_y_symmetry  - Whether to use the symmetry in y direction and only get indices for one slit.
%   polarisation    - Direction of polarisation. Use 'x', 'y' or 'z'.
%
% Outputs:
%   idx             - Canonical indices of slits at the left boundary (edges in polarisation direction).

    y_slit = [(-d - delta)/2, (-d + delta)/2, (d - delta)/2, (d + delta)/2]; % y values of upper and lower slit.
    for i = 1:length(y_slit)
        % Find y-index closest to actual y_slit value
        [~,y_idx(i)] = min(abs(msh.ymesh - y_slit(i)));
    end

    % Find all y-indices between slits
    if use_y_symmetry
        y_idx = y_idx(3):y_idx(4);
    else
        y_idx = [y_idx(1):y_idx(2), y_idx(3):y_idx(4)];
    end

    % Transform y indices to canonical indices
    if polarisation == 'x'
        idx = 1 + msh.nx * (y_idx-1);
    elseif polarisation == 'y'
        idx = 1 + msh.nx * (y_idx-1) + msh.np;
    elseif polarisation == 'z'
        idx = 1 + msh.nx * (y_idx-1) + 2*msh.np;
end