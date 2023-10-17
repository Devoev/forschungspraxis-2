function [e_pos, h_pos, s_pos] = analytic_sol_ypol(x_pos, E1, E2, lambda1, lambda2, d, a, n, daY)
%ANALYTIC_SOL calculates the analytic solution for the plane wave intensity
% at selected x-values (for y-polarized waves) 
% for a thin film interference of two y-polarised plane waves.
%INPUT:
%   x_pos - vector of x Positions on which the intensity is calculated
%   E1 - Field strength of first excitation
%   E2 - Field strnegth of second excitation
%   lambda1 - wavelength of first excitation
%   lambda2 - wavelength of second excitation
%   d - distance from excitation to thin film border
%   a - thickness of the thin film
%   n - refraction index of the thin film
%   daY - size of y surfaces in the model domain
%OUTPUT:
%   s_pos - analytic evaluations of poyntin vec in y direction 
%           => FIELD INTENSITY at x_pos
%                     a  
%                   |---|     
%   -->             * * *
%   -->             * * *
%   -->         1   * 2 *   3
%   -->             * * *
%   -->             * * *
%      |------------|
%     x=0          x=d       

dim = length(x_pos);
% init return values
e_pos = zeros(dim,1);
h_pos = zeros(dim,1);
s_pos = zeros(dim,1);

c = 3e8;
eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
% wave impedance for each area
Z1 = sqrt(mu0/eps0);
Z2 = sqrt(mu0/(n^2*eps0));
Z3 = Z1;
% wave numbers for first excitation
k1_l1 = (2*pi)/lambda1;
k2_l1 = (2*pi)/(lambda1/n);
k3_l1 = k1_l1;
% wave numbers for second excitation
k1_l2 = (2*pi)/lambda2;
k2_l2 = (2*pi)/(lambda2/n);
k3_l2 = k1_l2;

%% calculate coefficients for both excitations
[R1_l1, T2_l1, R2_l1, T3_l1] = get_coefficients(k2_l1, k3_l1, a, Z1, Z2, Z3);
[R1_l2, T2_l2, R2_l2, T3_l2] = get_coefficients(k2_l2, k3_l2, a, Z1, Z2, Z3);

%% calculate analytic electric field wave solutions (y-Pol)
E1_l1 = @(x) E1*(exp(-1j*k1_l1*x) + R1_l1*exp(1j*k1_l1*x));
E2_l1 = @(x) E1*(T2_l1*exp(-1j*k2_l1*x) + R2_l1*exp(1j*k2_l1*x));
E3_l1 = @(x) E1*T3_l1*exp(-1j*k3_l1*x);

E1_l2 = @(x) E2*(exp(-1j*k1_l2*x) + R1_l2*exp(1j*k1_l2*x));
E2_l2 = @(x) E2*(T2_l2*exp(-1j*k2_l2*x) + R2_l2*exp(1j*k2_l2*x));
E3_l2 = @(x) E2*T3_l2*exp(-1j*k3_l2*x);

%% calculate analytic magnetic field wave solutions (z-Pol)
H1_l1 = @(x) (E1/Z1)*(exp(-1j*k1_l1*x) - R1_l1*exp(1j*k1_l1*x));
H2_l1 = @(x) (E1/Z2)*(T2_l1*exp(-1j*k2_l1*x) - R2_l1*exp(1j*k2_l1*x));
H3_l1 = @(x) (E1/Z3)*T3_l1*exp(-1j*k3_l1*x);

H1_l2 = @(x) (E2/Z1)*(exp(-1j*k1_l2*x) - R1_l2*exp(1j*k1_l2*x));
H2_l2 = @(x) (E2/Z2)*(T2_l2*exp(-1j*k2_l2*x) - R2_l2*exp(1j*k2_l2*x));
H3_l2 = @(x) (E2/Z3)*T3_l2*exp(-1j*k3_l2*x);

%% calculate intensity
% shift analytic solution to match model domain
eval_pos = x_pos - d;

for idx = 1:dim
    x_val = eval_pos(idx);
    if x_val<0
        % in area 1 -> before thin film
        e_pos(idx) = E1_l1(x_val) + E1_l2(x_val);
        h_pos(idx) = H1_l1(x_val) + H1_l2(x_val);
    elseif 0<=x_val && x_val<a
        % in area 2 -> thin film
        e_pos(idx) = E2_l1(x_val) + E2_l2(x_val);
        h_pos(idx) = H2_l1(x_val) + H2_l2(x_val);
    else
        % in area 3 -> after thin film
        e_pos(idx) = E3_l1(x_val) + E3_l2(x_val);
        h_pos(idx) = H3_l1(x_val) + H3_l2(x_val);
    end
    % poynting vector in x direction:
    s_pos(idx) = (1/2)*e_pos(idx)*conj(h_pos(idx));
end

    %% Solving the plane wave LG
    function [R1, T2, R2, T3] = get_coefficients(k2, k3, a, Z1, Z2, Z3)
        % entry impedance
        Ze = Z3*(1+Z2/Z3*tanh(1j*k2*a))/(1+Z3/Z2*tanh(1j*k2*a));
        % calc the reflection and transmission coefficients
        R1 = (Ze-Z1)/(Ze+Z1);

        T2 = Z2*(R1-1)/(2*Z1) + 1/2 + R1/2;
        R2 = 1+R1-T2;
        
        T3 = (T2*exp(-1j*k2*a)+R2*exp(1j*k2*a))*exp(1j*k3*a);
    end
end