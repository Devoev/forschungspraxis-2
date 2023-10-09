function I = calc_intensity_fd(ebow_abs)
% CALC_INTENSITY Calculates the intensity of a complex ebow in FD.
%
% Inputs:
%   ebow_abs    - Absolute value of integrated electrical field. Vector of size (np) or smaller.
%
% Outputs:
%   I       - Intensity. Vector of size(ebow_abs).

    % Constants
    c = 3e8;
    eps = 8.854e-12;

    % Intensity calculation
    I = c*eps/4 * ebow_abs.^2;              % Intensity formula
end