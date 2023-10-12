
function [I, bright, dark] = intensityCalcSignleSlit(Imax, d, l, lambda, y)

    % Calculate intensity over y-coordinates
    I = Imax * (sinc(d/lambda * sin(atan(y/l)))).^2;
    
    % Calculate maximum theta in the displayed area
    theta_max = atan(max(y)/l);
    
    % Calculate maximum number of dark fringes
    m = floor(sin(theta_max) * d / lambda);
    
    % Calculate y-values of dark fringes
    dark = tan(asin((1:m)*lambda/d))*l;
    
    % Calculate y-values of bright fringes
    bright = dark(1:end-1) + diff(dark)*1/2;

end
