function I = helmholtz_ana(E0, lambda, d, delta, L, y, n)
% HELMHOLTZ_ANA Calculates the intensity at the screen using the analytical solution of the Helmholtz equation.
%
% Inputs:
%   E0      - Max field strength.
%   lambda  - Wavelength.
%   d       - Slit distance.
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
    y_slit_1 = linspace((-d - delta)/2, (-d + delta)/2, n);
    y_slit_2 = linspace((d - delta)/2, (d + delta)/2, n);
    for i = 1:n
        y1 = y - y_slit_1(i);
        y2 = y - y_slit_2(i);
        r1 = sqrt(L^2 + y1.^2);
        r2 = sqrt(L^2 + y2.^2);
        E_tot = E_tot + E(r1) + E(r2);
    end

    I = c*eps/2 * abs(E_tot).^2; % Add factor of 1/2 ?
end