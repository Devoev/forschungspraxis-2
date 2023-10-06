function [e_pos] = analytic_sol(x_pos, E1, E2, lambda1, lambda2, d, a, n)
%ANALYTIC_SOL calculates the analytic solution at selected x-values 
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
%OUTPUT:
%   e_pos - analytic evaluations of el. field at all x positions (x_pos)
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
e_pos = zeros(dim,1);

eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
% wave impedance for each area
Z1 = sqrt(mu0/eps0);
Z2 = sqrt(mu0/n^2*eps0);
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

%% calculate analytic wave solutions
E1_l1 = @(x) E1*(exp(-1j*k1_l1*x) + R1_l1*exp(1j*k1_l1*x));
E2_l1 = @(x) E1*(T2_l1*exp(-1j*k2_l1*x) + R2_l1*exp(1j*k2_l1*x));
E3_l1 = @(x) E1*T3_l1*exp(-1j*k3_l1*x);

E1_l2 = @(x) E2*(exp(-1j*k1_l2*x) + R1_l2*exp(1j*k1_l2*x));
E2_l2 = @(x) E2*(T2_l2*exp(-1j*k2_l2*x) + R2_l2*exp(1j*k2_l2*x));
E3_l2 = @(x) E2*T3_l2*exp(-1j*k3_l2*x);

%% loop over x positions and calculate 
% correct analytic x offset
eval_pos = x_pos - d;
for idx = 1:dim
    x_val = eval_pos(idx);
    if x_val<0
        % in layer 1
        e_pos(idx) = E1_l1(x_val) + E1_l2(x_val);
    elseif 0<=x_val && x_val<a
        % in layer 2
        e_pos(idx) = E2_l1(x_val) + E2_l2(x_val);
    else
        % in layer 3
        e_pos(idx) = E3_l1(x_val) + E3_l2(x_val);
    end
end

    function [R1, T2, R2, T3] = get_coefficients(k2, k3, a, Z1, Z2, Z3)
        % helper term
        A = 2*exp(-1j*k2*a)*(Z3-Z2)/(Z3+Z2);
        % ToDo: calc all reflection and transmission coefficients
        R2 = (2*A*Z2)/((Z1+Z2)*(1-(Z1-Z2)*A/(Z1+Z2)));
        T2 = (2*Z2 + R2*(Z1-Z2))/(Z1+Z2);
        R1 = T2+R2-1;
        T3 = (T2*exp(-1j*k2*a)+R2*exp(1j*k2*a))*exp(1j*k3*a);
    end

end