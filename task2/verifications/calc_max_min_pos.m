function [d_max,d_min] = calc_max_min_pos(y,L,d,lambda)
% calc_max_min_pos Calculates the position of minima and maxima on the screen.
%
% Inputs:
%   y         - y values for the screen
%   L         - Length of Domain
%   d         - Distance of the 2 slits
%   lambda    - Wavelength
%
% Outputs:
%   d_max           - position of maxima on screen
%   d_min           - position of minima on screen

    j=1;
    for k = ceil(min(y)/(lambda*(L/d))):floor(max(y)/(lambda*(L/d)))
        d_min(j) = (k*lambda)*(L/d);
        d_max(j) = ((2*k+1)*lambda/2)*(L/d);
        j = j+1;
    end
end