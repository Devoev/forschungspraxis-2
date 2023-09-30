function I = intensity_farfield(E0, lambda, delta, L, y)
% INTENSITY_FARFIELD Calculates the intensity at the screen using the analytical double slit formula.
% Adapted for 1 slit
%
% Inputs:
%   E0      - Max field strength.
%   lambda  - Wavelength.
%   d       - Slit distance.
%   delta   - Slit width.
%   L       - Screen distance.
%   y       - Array of y coordinates on the screen.
%
% Outputs:
%   I       - Intensity. Array of length(y).

    % Constants
    c = 3e8;
    eps = 8.854e-12;

    % Variables
    theta = atan(y/L);
    I0 = c*eps*2*E0^2;


    % w/ diffraction
    %x = pi*delta*sin(theta)/lambda;
    I = I0 *( (sin(pi*delta*sin(theta) / lambda) ./ (pi*delta*sin(theta) / lambda) )).^2;
end