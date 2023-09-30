function [Meps_s, Mmui_s] = calcpml_2D(msh, bc, Meps, Mmui)
% calcpml2d. von A.G.
% 
% Inputs:
%   NGRID   - mesh dimentions list with [nx, ny]
%   NPML    - number of mesh layers [L1, L2, L3, L4]
%
%          L3 (top)               y
%       +++++++               o------->
%   L4  +     +  L2           |
% (src) +     + (screen)    x | 
%       +++++++               |
%          L1 (btm)           v
%
% Outputs:
%   Meps_s   - permeability matrix with PML boundary material
%   Mmui_s    - inverse permittivity matrix with PML boundary material

% pml settings
sigma_max = 1;
a_max = 3;
p = 3;

% input
% TODO - unify coordinates x:=x y:=y
Nx = msh.ny;
Ny = msh.nx;

L1 = bc.NPML(2); 
L2 = bc.NPML(1);
L3 = bc.NPML(4);
L4 = bc.NPML(3);

% Calculate vacuum wave impedance
if any(ismember(fieldnames(bc),"epsilon0"))
    imp0 = sqrt(1/bc.mu0i/bc.epsilon0);
else
    imp0 = 376.73;
end

% init tensors ------------------------------------------------------------
sx = sparse(ones(Nx, Ny));
sy = sparse(ones(Nx, Ny));

%% calc sx tensor comp ----------------------------------------------------
% add L2 direction PML
if L2 > 1
    for nx = 1:L2
        v = nx/L2;
        % set column from inside to outside
        sx(L2-nx+1,:) = impedance(v, L2, sigma_max, a_max, p, imp0); 
    end
end
% add L4 direction PML
if L4 > 1
    for nx = 1:L4
        v = nx/L4;
        % set column from inside to outside
        sx(Nx-L4+nx,:) = impedance(v, L4, sigma_max, a_max, p, imp0); 
    end
end
%% calc sy tensor comp ----------------------------------------------------
% add L1 direction PML
if L1 > 1
    for ny = 1:L1
        v = ny/L1;
        % set row from inside to outside of boundary
        sy(:,Ny-L1+ny) = impedance(v, L1, sigma_max, a_max, p, imp0);
    end
end
% add +Y direction PML
if L3 > 1
    for ny = 1:L3
        v = ny/L3;
        % set row from inside to outside of boundary
        sy(:,L3-ny+1) = impedance(v, L3, sigma_max, a_max, p, imp0);
    end
end

%% build and apply the 3np x 3np S tensor

sx_v = reshape(sx', [], 1);
sy_v = reshape(sy', [], 1);
S = sparse(diag([sy_v./sx_v; sx_v./sy_v; sx_v.*sy_v]));

Meps_s = S*Meps;
Mmui_s = S*Mmui;
    
% calculate matching impedance for a certain pml layer --------------------
function val = impedance(v, Lv, sigma_max, a_max, p, imp0)

    % calc sigma
    sigma = sigma_max*sin(pi*v/(2*Lv))^2;
    % calc attenuation
    a = 1+a_max*(v/Lv)^p;
    
    val = a*(1+1j*imp0*sigma);
    
end

end