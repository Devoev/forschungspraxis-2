function I = calc_intensity_td(ebow)
% CALC_INTENSITY_TD Calculates the time averaged intensity of ebow in TD.
%
% Inputs:
%   ebow    - Integrated electrical field over time. Matrix of size (3*np,nt).
%
% Outputs:
%   I       - Intensity. Vector of size(3*np).

    % Constants
    c = 3e8;
    eps = 8.854e-12;

    [n, nt] = size(ebow);

    % Intensity formula
    I = c*eps/4 * ebow(0).^2;   % TODO: Add time averaging
end