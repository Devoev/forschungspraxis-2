function [MAT] = conductivePML_2D(bc, msh, MAT, f)
% conductivePML_2D creates a conductivity matrix, creating an absorbing
% boundary for a given domain
%
% Input
% bc                -Boundary condition object, including the information
%                    where to set an open boundary and the size of the PML
%                    layer in cells. Example:
%                    bc.bc = ["OPEN", "OPEN", "OPEN", "OPEN"];
%                    bc.NPML = [1,1,1,1]*6*elem_per_wavelength;
%                    Recommended is a width of 6 wavelengths per PML layer
%                    regarding the longest wavelength simulated.
% msh               -2D mesh object
% MAT               -MAT object containg all relevant matrices
% f                 -Excitation frequency
%
% Output
% MAT               -MAT object, now including a conductivity matrix
%                    to simulate the absorbing boundary


% Get basic mesh parameters
nx = msh.nx;
ny = msh.ny;

% Calculate basic value for PML conductivity
kappa = 0.01 * 2 * pi * f * MAT.epsilon0;

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
