function I = intensity_ana(E0, lambda, d, L, y)
% INTENSITY_ANA Calculates the intensity at the screen using the analytical double slit formula.
%
% Inputs:
%   E0      - Max field strength.
%   lambda  - Wavelength.
%   d       - Slit distance.
%   L       - Screen distance.
%   y       - Array of y coordinates on the screen.
%
% Outputs:
%   I       - Intensity. Array of length(y).

    I = 4*E0^2 * cos(pi*d/(lambda*L)*y).^2;
end