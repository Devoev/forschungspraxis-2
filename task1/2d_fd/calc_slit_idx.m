function idx = calc_slit_idx(msh, d, delta)
% CALC_SLIT_IDX Calculates the indices for the the 2 slits of the double slit.
%
% Inputs:
%   msh     - Mesh.
%   d       - Slit distance.
%   delta   - Slit width.
%
% Outputs:
%   idx     - Canonical indices of slits at the left boundary.

    y_slit = [(-d - delta)/2, (-d + delta)/2, (d - delta)/2, (d + delta)/2]; % y values of upper and lower slit.
    for i = 1:length(y_slit)
        % Find y-index closest to actual y_slit value
        [~,y_idx(i)] = min(abs(msh.ymesh - y_slit(i)));
    end
    y_idx = [y_idx(1):y_idx(2), y_idx(3):y_idx(4)]; % Find all y-indices between slits
    idx = msh.nx * (y_idx-1); % Transform y indices to canonical indices
end