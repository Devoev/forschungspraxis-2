function [intensity] = calc_intensity(x_pos, E1, E2, lambda1, lambda2, L, a, n)
%CALC_INTENSITY calculates the analytic values for the resulting wave
%intensity at a list of x-positions.
%INPUT:
%   x_pos - vector of x Positions on which the intensity is calculated
%   E1 - Field strength of first excitation
%   E2 - Field strnegth of second excitation
%   lambda1 - wavelength of first excitation
%   lambda2 - wavelength of second excitation
%   L - distance from excitation to thin film border
%   a - thickness of the thin film
%   n - refraction index of the thin film layer
%OUTPUT:
%   intensity - vector of intensity values of same length as x_pos input
%                     a  
%                   |---|     
%   -->             * * *
%   -->             * * *
%   -->         0   * 1 *   2
%   -->             * * *
%   -->             * * *
%      |------------|-----------|
%     x=0         x=L/2        x=L

eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
% for first excitation
z0 = sqrt(mu0/eps0);
k0 = (2*pi)/lambda1;
z1 = sqrt(mu0/(eps0*n^2));
k1 = (2*pi)/(lambda1/n);
z2 = z0;
k2 = k0;

% -> border 0 to 1 coefficients
r0 = (z1-z0)/(z1+z0);
d0 = 1;
% -> border 1 to 2





end