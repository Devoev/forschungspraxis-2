function [I,y] = calc_intensity(msh, ebow, offset)
% CALC_INTENSITY Calculates the intensity at the screen using the numerical result ebow.
%
% Inputs:
%   msh     - Mesh struct.
%   ebow    - Phasor of integrated electrical field. Vector of size (3*np).
%   offset  - Offsets from the boundary in each direction. Vector of size (4).
%
% Outputs:
%   I       - Intensity. Vector of size(ny - offset(2) - offset(4)).
%   y       - Y coordinates at which the intensity is evaluated. Vector of size (I).

    % Y coordinates
    idx_yoffset = 1+offset(1):length(msh.ymesh)-offset(3);  % y indices at the screen with offset
    y = msh.ymesh(idx_yoffset);

    % Intensity calculation
    idx_screen = msh.nx * (1:msh.ny) - offset(2);           % Canonical indices at the screen
    ebow_abs = calc_abs_field(msh, ebow);                   % Calculate abs value
    ebow_abs = ebow_abs(idx_screen);                        % Evaluate field at screen
    ebow_abs = ebow_abs(idx_yoffset);                       % Offset from y boundary
    I = calc_intensity_fd(ebow_abs);                        % Intensity formula
end