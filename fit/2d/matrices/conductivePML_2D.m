function [MAT] = conductivePML_2D(bc, msh, MAT, fmin)
% createMeps_2D creates permittivity matrix for 2D mesh
%
% Input
% msh               -2D mesh object
% ds                -primary edge matrix
% da                -primary face matrix
% dat               -dual face matrix
% rel_eps           -vector with distribution of relative permeability
% eps0              -permeability of vacuum
%
% Output
% meps              -permittivity matrix


% Get basic mesh parameters
nx = msh.nx;
ny = msh.ny;

% Calculate basic value for PML conductivity
kappa = 0.01 * 2 * pi * fmin * MAT.epsilon0;

% Counter for kth box
k = 1;

% Initiate conductivity for associated PML-layer for boundary at y = ymin
if bc.open_bc(1)
    for i = 1:1:bc.NPML(1)
        boxesKappa(k).box = [1, nx, bc.NPML(1)+1-i, bc.NPML(1)+2-i];  %#ok<*AGROW> 
        boxesKappa(k).value = (i/bc.NPML(1)*20 * kappa);  
        k = k + 1;
    end
end

% Initiate conductivity for associated PML-layer for boundary at x = xmax
if bc.open_bc(2)
    for i = 1:1:bc.NPML(2)
        boxesKappa(k).box = [nx-bc.NPML(2)-1+i, nx-bc.NPML(2)+i, 1, ny];  
        boxesKappa(k).value = (i/bc.NPML(2)*20 * kappa);  
        k = k + 1;
    end
end

% Initiate conductivity for associated PML-layer for boundary at y = ymax
if bc.open_bc(3)
    for i = 1:1:bc.NPML(3)
        boxesKappa(k).box = [1, nx, ny-bc.NPML(3)-1+i, ny-bc.NPML(3)+i];  
        boxesKappa(k).value = (i/bc.NPML(3)*20 * kappa);  
        k = k + 1;
    end
end

% Initiate conductivity for associated PML-layer for boundary at x = xmin
if bc.open_bc(4)
    for i = 1:1:bc.NPML(4)
        boxesKappa(k).box = [bc.NPML(4)+1-i, bc.NPML(4)+2-i, 1, ny];  
        boxesKappa(k).value = (i/bc.NPML(4)*20 * kappa);  
        k = k + 1;
    end
end


kap = boxMesher_2D(msh, boxesKappa, 0);
MAT.kaps = createMeps_2D(msh, MAT.ds, MAT.da, MAT.dat, kap, 1);

end
