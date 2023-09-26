function res = fresnel_number(delta, L, lambda)
%% FRESNEL_NUMBER Calculates the fresnel number for the given parameters.
%
% Inputs:
%   delta   - Slit with.
%   L       - Screen distance.
%   lambda  - Wavelength.
%
% Outputs:
%   res     - Value of the fresnel number.

    res = delta^2/(L*lambda);
    assert(res < .2, "Value of fresnel number is " + res);
end