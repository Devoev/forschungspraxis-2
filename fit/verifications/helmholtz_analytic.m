function [E_ana] =helmholtz_analytic(lambda, x, y, elem_per_wavelength, omega, excitation)

kappa = (2*pi)/lambda; %Wavenumber

%screen: 1 lambda before end of domain
%d = distance of elements of the screen to source
d_x = (x-1)*lambda;
d_y = linspace(-(0.5*y-1)*lambda,(0.5*y-1)*lambda,(elem_per_wavelength*(x-2)+1));
d = sqrt(d_x .* d_x + d_y .* d_y);

E_komplex =(excitation*exp(-1i*kappa*d) ./d);
%E_real = real(E_komplex*exp(-1i*omega));
E_ana =real(E_komplex).^2;

end