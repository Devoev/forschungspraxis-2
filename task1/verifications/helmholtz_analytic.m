function [E_ana] =helmholtz_analytic(lambda, x, ymesh, elem_per_wavelength, excitation, yidx, h)

% Inputs:
% lambda: Wavelength; x: Distance of screen in x dir; ymesh: mesh in y dir.
% elem_per_wavelength: elements per wavelength, excitation: Amplitude of
% excitation, yidx: Y-Index of Sources, h: clac domain in y dir.

kappa = (2*pi)/lambda; % Wavenumber

% d = distance of elements of the screen to source
d_x = x;

E_sum = zeros(1,length(ymesh));

for i =1 : length(yidx)     % go through all excitations
    d_y =ymesh - (-h/2 + (yidx(i)/(h/lambda*elem_per_wavelength)) * h);
    d = sqrt(d_x .* d_x + d_y .* d_y);

    E_komplex =(excitation*exp(-1i*kappa*d) ./d);
    E_sum = E_sum + E_komplex;
    %plot(ymesh, real(E_sum))   % reference plot for the last excitation
end
E_ana =real(E_sum).^2;
end