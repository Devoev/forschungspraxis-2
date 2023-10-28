function [S] = CalcPoyntinvectorXY(msh, ebow, hbow, ds, dst)
% CalcPowerSurfaceXY calculates the power through each surface in x- and
% y-direction utilizing the Poyntinvector. Applyable for time and frequency
% domain.
%
% Inputs:
%   msh                 - Mesh struct
%   ebow                - Vector with integrated electric field
%   hbow                - Vector with integrated magnetic field
%   ds                  - Primary edge matrix
%   dst                 - Dual edge matrix
%
% Outputs:
%   S                   - Integrated Poyntinvector with entries for faces
%                         x- and y-direction. Equals the power through each
%                         surface


%% Get basic mesh parameters
nx = msh.nx;
ny = msh.ny;
np = msh.np;
Mx = msh.Mx;
My = msh.My;


%% Get electric and magnetic field quantities
E = nullInv(ds) * ebow;
H = nullInv(dst) * hbow;


%% Calculate power through each relevant surface in y-direction (common)

% Indices for common calculation (majority of the mesh)
n_subdom = 1 + ((1:nx-1)-1) * Mx;
n_subdom = repmat(n_subdom, ny-2, 1)' + My * (1:ny-2);
n_subdom = reshape(n_subdom, (nx-1)*(ny-2), 1);

% Calculate z-component of the electric field at each surface
Ez = 1/2 * (E(n_subdom+2*np) + E(n_subdom+2*np+Mx));

% Calculate x-component of the magnetic field at each surface
Hx = 1/4 * (H(n_subdom) + H(n_subdom+Mx) + H(n_subdom-My) + H(n_subdom-My+Mx));

% Calculate x-component of the electric field at each surface
Ex = E(n_subdom);

% Calculate z-component of the magnetic field at each surface
Hz = 1/2 * (H(n_subdom+2*np) + H(n_subdom+2*np-My));

% Calculate Poyntinvector to each surface
Sy = zeros(np,1);
Sy(n_subdom) = Ez .* conj(Hx) - Ex .* conj(Hz);


%% Calculate power through each relevant surface in y-direction (y = ymin)

% Indices for surfaces on boundary y = ymin
n_subdom = 1 + ((1:nx-1)-1) * Mx;

% Calculate z-component of the electric field at each surface
Ez = 1/2 * (E(n_subdom+2*np) + E(n_subdom+2*np+Mx));

% Calculate x-component of the magnetic field at each surface
Hx = 1/2 * (H(n_subdom) + H(n_subdom+Mx));

% Calculate x-component of the electric field at each surface
Ex = E(n_subdom);

% Calculate z-component of the magnetic field at each surface
Hz = H(n_subdom+2*np);

% Calculate Poyntinvector to each surface
Sy(n_subdom) = Ez .* conj(Hx) - Ex .* conj(Hz);


%% Calculate power through each relevant surface in y-direction (y = ymax)

% Indices for surfaces on boundary y = ymax
n_subdom = 1 + ((1:nx-1)-1) * Mx + (ny-1) * My;

% Calculate z-component of the electric field at each surface
Ez = 1/2 * (E(n_subdom+2*np) + E(n_subdom+2*np+Mx));

% Calculate x-component of the magnetic field at each surface
Hx = 1/2 * (H(n_subdom-My) + H(n_subdom+Mx-My));

% Calculate x-component of the electric field at each surface
Ex = E(n_subdom);

% Calculate z-component of the magnetic field at each surface
Hz = H(n_subdom+2*np-My);

% Calculate Poyntinvector to each surface
Sy(n_subdom) = Ez .* conj(Hx) - Ex .* conj(Hz);


%% Calculate power through each relevant surface in x-direction (common)

% Indices for common calculation (majority of the mesh)
n_subdom = 1 + ((2:nx-1)-1) * Mx;
n_subdom = repmat(n_subdom, ny-1, 1)' + My * (0:ny-2);
n_subdom = reshape(n_subdom, (nx-2)*(ny-1), 1);

% Calculate y-component of the electric field at each surface
Ey = E(n_subdom+np);

% Calculate z-component of the magnetic field at each surface
Hz = 1/2 * (H(n_subdom+2*np) + H(n_subdom+2*np-Mx));

% Calculate z-component of the electric field at each surface
Ez = 1/2 * (E(n_subdom+2*np) + E(n_subdom+2*np+My));

% Calculate y-component of the magnetic field at each surface
Hy = 1/4 * (H(n_subdom+np) + H(n_subdom+np-Mx) + H(n_subdom+np+My) + H(n_subdom+np+My-Mx));

% Calculate Poyntinvector to each surface
Sx = zeros(np,1);
Sx(n_subdom) = Ey .* conj(Hz) - Ez .* conj(Hy);


%% Calculate power through each relevant surface in x-direction (x = xmin)

% Indices for surfaces on boundary x = xmin
n_subdom = 1 + ((1:ny-1)-1) * My;

% Calculate y-component of the electric field at each surface
Ey = E(n_subdom+np);

% Calculate z-component of the magnetic field at each surface
Hz = H(n_subdom+2*np);

% Calculate z-component of the electric field at each surface
Ez = 1/2 * (E(n_subdom+2*np) + E(n_subdom+2*np+My));

% Calculate y-component of the magnetic field at each surface
Hy = 1/2 * (H(n_subdom+np) + H(n_subdom+np+My));

% Calculate Poyntinvector to each surface
Sx(n_subdom) = Ey .* conj(Hz) - Ez .* conj(Hy);


%% Calculate power through each relevant surface in x-direction (x = xmax)

% Indices for surfaces on boundary x = xmax
n_subdom = 1 + (nx-1) * Mx + ((1:ny-1)-1) * My;

% Calculate y-component of the electric field at each surface
Ey = E(n_subdom+np);

% Calculate z-component of the magnetic field at each surface
Hz = H(n_subdom+2*np-Mx);

% Calculate z-component of the electric field at each surface
Ez = 1/2 * (E(n_subdom+2*np) + E(n_subdom+2*np+My));

% Calculate y-component of the magnetic field at each surface
Hy = 1/2 * (H(n_subdom+np-Mx) + H(n_subdom+np+My-Mx));

% Calculate Poyntinvector to each surface
Sx(n_subdom) = Ey .* conj(Hz) - Ez .* conj(Hy);


%% Set up complete Poyntinvector and apply TD/FD

% Set up result vector
S = [Sy; Sx; zeros(np,1)];

% Multiply with 0.5, if imaginary
if isreal(S) == false
    S = 1/2 *S;
end
