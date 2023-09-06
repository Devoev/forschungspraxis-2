function res = fresnel_number(delta, h, lambda)
%% FRESNEL_NUMBER Calculates the fresnel number for the given parameters.
%
% Inputs:
%   delta   - Slit with.
%   h       - Screen distance.
%   lambda  - Wavelength.
%
% Outputs:
%   res     - Value of the fresnel number.

    res = delta^2/(4*h*lambda);
    assert(res < .2);
end