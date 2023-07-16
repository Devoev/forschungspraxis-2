function [sx, sy] = calcpml2D(NGRID, NPML)
% calcpml2d. von A.G.
% 
% Inputs:
%   NGRID   - mesh dimentions list with [nx, ny]
%   NPML    - number of mesh layers [L1, L2, L3, L4]
%
%          L2 (top)               y
%       +++++++               o------->
%   L3  +     +  L1           |
% (src) +     + (screen)    x | 
%       +++++++               |
%          L4 (btm)           v
%
% Outputs:
%   sx    - x pml tensor w. dimension like mesh Nx x Ny
%   sy    - y pml tensor w. dimension like mesh Nx x Ny

% pml settings
sigma_max = 1;
a_max = 3;
p = 3;

% input
% TODO - unify coordinates x:=x y:=y
Nx = NGRID(2);
Ny = NGRID(1);

L1 = NPML(2); 
L2 = NPML(3);
L3 = NPML(4);
L4 = NPML(1);

% init tensors ------------------------------------------------------------
sx = sparse(ones(Nx, Ny));
sy = sparse(ones(Nx, Ny));

% calc sx tensor ----------------------------------------------------------
% add L2 direction PML
if L2 > 1
    for nx = 1:L2
        v = nx/L2;
        % set column from inside to outside
        sx(L2-nx+1,:) = impedance(v, L2, sigma_max, a_max, p); 
    end
end
% add L4 direction PML
if L4 > 1
    for nx = 1:L4
        v = nx/L4;
        % set column from inside to outside
        sx(Nx-L4+nx,:) = impedance(v, L4, sigma_max, a_max, p); 
    end
end
% calc sy tensor ----------------------------------------------------------
% add L1 direction PML
if L1 > 1
    for ny = 1:L1
        v = ny/L1;
        % set row from inside to outside of boundary
        sy(:,Ny-L1+ny) = impedance(v, L1, sigma_max, a_max, p);
    end
end
% add +Y direction PML
if L3 > 1
    for ny = 1:L3
        v = ny/L3;
        % set row from inside to outside of boundary
        sy(:,L3-ny+1) = impedance(v, L3, sigma_max, a_max, p);
    end
end

    
% calculate matching impedance for a certain pml layer --------------------
function val = impedance(v, Lv, sigma_max, a_max, p)
    imp0 = 376.73;
    % calc sigma
    sigma = sigma_max*sin(pi*v/(2*Lv))^2;
    % calc attenuation
    a = 1+a_max*(v/Lv)^p;
    
    val = a*(1+1j*imp0*sigma);
end

end