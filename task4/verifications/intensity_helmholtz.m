function I = intensity_helmholtz(E0, lambda, delta, L, y, n)
% INTENSITY_HELMHOLTZ Calculates the intensity at the screen using the analytical solution of the Helmholtz equation.
% Adapted for 1 slit.
%
% Inputs:
%   E0      - Max field strength.
%   lambda  - Wavelength.
%   delta   - Slit width.
%   L       - Screen distance.
%   y       - Array of y coordinates on the screen.
%   n       - Number of discrete excitations in each slit.
%
% Outputs:
%   I       - Intensity. Array of length(y).

    % Constants
    c = 3e8;
    eps = 8.854e-12;
    kappa = 2*pi/lambda; % Wavenumber

    E_tot = zeros(1,length(y));
    E = @(r) E0*exp(-1i*kappa*r) ./ (kappa*r);

    % Iteration over all excitations
    y_slit_1 = linspace(-delta/2, delta/2, n);
    for i = 1:n
        y1 = y - y_slit_1(i);
        r1 = sqrt(L^2 + y1.^2);
        E_tot = E_tot + E(r1);
    end

    I = c*eps/2 * abs(E_tot).^2; % Add factor of 1/2 ?
end