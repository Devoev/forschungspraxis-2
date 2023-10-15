function [E_ana, H_ana] = IncidentPlaneWaveAnalytic(positions, E0, f, epsr1, epsr2, MAT)

% Find index position = 0 for medium change
idx_x0 = find(positions == 0);

% Offset for analytic solution
positions = positions - positions(1);

% Calculate wave impedances for both media
Z1 = sqrt(1/MAT.mu0i/MAT.epsilon0/epsr1);
Z2 = sqrt(1/MAT.mu0i/MAT.epsilon0/epsr2);

% Calculate transmission and refelction coefficients 
t = 2 * Z2 / (Z1 + Z2);
r = (Z2 - Z1) / (Z2 + Z1);

% Calculate propagation constants for both media
k1 = 2 * pi * f * sqrt(MAT.epsilon0 * epsr1 / MAT.mu0i);
k2 = 2 * pi * f * sqrt(MAT.epsilon0 * epsr2 / MAT.mu0i);

% Calculate electric and magnetic field for media 1
E_ana_1 = E0 * exp(-1j * k1 * positions) - r * E0 * exp(1j * k1 * positions);
H_ana_1 = E0/Z1 * exp(-1j * k1 * positions) + r * E0/Z1 * exp(1j * k1 * positions);

% Calculate electric and magnetic field for media 2
E_ana_2 = -t * E0 * exp(-1j * k2 * positions);
H_ana_2 = -t * E0/Z2 * exp(-1j * k2 * positions);

% Combine solutions
E_ana = E_ana_1;
E_ana(idx_x0:end) = E_ana_2(idx_x0:end);
H_ana = H_ana_1;
H_ana(idx_x0:end) = H_ana_2(idx_x0:end);

end
