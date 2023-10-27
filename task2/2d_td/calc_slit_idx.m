function idx = calc_slit_idx(msh, d, delta, polarisation)
% CALC_SLIT_IDX Calculates the indices for the the slit of task 2.
%
% Inputs:
%   msh             - Mesh.
%   d               - Slit distance.
%   delta           - Slit width.
%   polarisation    - Direction of polarisation. Use 'x', 'y' or 'z'.
%
% Outputs:
%   idx             - Canonical indices of slits at the left boundary (edges in polarisation direction).

    y_slit = [(d - delta)/2, (d + delta)/2]; % y values of upper and lower slit.
    for i = 1:length(y_slit)
        % Find y-index closest to actual y_slit value
        [~,y_idx(i)] = min(abs(msh.ymesh - y_slit(i)));
    end

    % Find all y-indices between slits
    y_idx = y_idx(1):y_idx(2);

    % Transform y indices to canonical indices
    if polarisation == 'x'
        idx = msh.nx * (y_idx-1);
    elseif polarisation == 'y'
        idx = msh.nx * (y_idx-1) + msh.np;
    elseif polarisation == 'z'
        idx = msh.nx * (y_idx-1) + 2*msh.np;
end