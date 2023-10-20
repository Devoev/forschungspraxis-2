function [S1, S3] = AnaSolPoyntin(E0, f, epsilon_1, epsilon_2, mui_1, mui_2, d)

% Calculate wave impedances of all spaces
Z1 = sqrt(1/epsilon_1/mui_1);
Z2 = sqrt(1/epsilon_2/mui_2);
Z3 = Z1;

% Calculate propagation constant inside the film
k2 = 2 * pi * f * sqrt(epsilon_2/mui_2);

% Calculate impedance for transition in third space
Ze = Z3 * (1 + Z2 / Z3 * 1j * tan(k2 * d)) / (1 + Z3 / Z2 * 1j * tan(k2 * d));

% Calculate magnitudes of reflection and transmission coefficients
r = abs((Ze - Z1) / (Ze + Z1));
d = abs(2 * Ze / (Ze + Z1));

% Calculate Poyntinvector in space before the film
S1 = 0.5 * E0 * conj(E0) / Z1 * (1 - r^2);

% Calculate Poyntinvector in space after the film
S3 = 0.5 * d^2 * E0 * conj(E0) / Z3;

end
