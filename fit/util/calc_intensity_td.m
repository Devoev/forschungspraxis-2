function I = calc_intensity_td(ebow_abs)
% CALC_INTENSITY_TD Calculates the time averaged intensity of ebow in TD.
%
% Inputs:
%   ebow_abs    - Absolute value of integrated electrical field over time. Matrix of size (n,nt) with some n<=np.
%
% Outputs:
%   I           - Intensity. Vector of size(n).

    [n, nt] = size(ebow_abs);
    ebow_abs = full(ebow_abs);

    % Calculate intensity over time
    It = zeros(n,nt);
    for i=1:nt
        It(:,i) = calc_intensity_fd(ebow_abs(:,i));
    end

    % Average intensity
    I = arrayfun(@(i) mean(It(i,:)), 1:n);
end