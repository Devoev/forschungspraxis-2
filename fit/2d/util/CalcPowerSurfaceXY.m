function [S] = CalcPowerSurfaceXY(msh, ebow, hbow, ds, dst, da)
% CalcPowerSurfaceXY calculates the power through each surface in x- and
% y-direction utilizing the Poyntinvector. Applyable for time and frequency
% domain
%
% Inputs:
%   msh                 - Mesh struct
%   bc                  - Boundary condition struct
%   material_regions    - Struct containing information about the material
%                         regions
%                         Example for regions for relative permittivity:
%                         boxesEpsilonR(1).box = [1, nx, 1, ny];
%                         boxesEpsilonR(1).value = 1;
%                         material_regions.boxesEpsilonR = boxesEpsilonR;
%                         -> Same ist needed for inverse permeability
%                         
%   operators           - List of desired operators like
%                         operators = ["CurlP", "SourceD"];
%                         -> Only operators named in the list are returned
%                         in MAT object. 
%                         CurlP -> Primary curl
%                         SourceP -> Primary source
%                         SourceD -> Dual source
%
% Outputs:
%   MAT                 - Struct containing all edited matrices


%% Get basic mesh parameters
nx = msh.nx;
ny = msh.ny;
np = msh.np;
Mx = msh.Mx;
My = msh.My;


%% Get electric and magnetic field quantities
E = nullInv(ds) * ebow;
H = nullInv(dst) * hbow;


%% Create extended vectors for calculation of the Pyntinvectors

% Widen actuel domain with one point in each direction
nx_ext = nx+2;
ny_ext = ny+2;
np_ext = nx_ext * ny_ext;
Mx_ext = 1;
My_ext = nx_ext;
E_ext = zeros(3*np_ext,1);
H_ext = zeros(3*np_ext,1);

% Calculate indices of entries of E and H in E_ext and H_ext
n_rel = 1 + ((2:nx_ext-1)-1)*Mx_ext + My_ext;
n_rel = reshape(repmat(n_rel,ny,1)' + (0:ny-1) * My_ext,np,1);
n = [n_rel;n_rel+np_ext;n_rel+2*np_ext];

% Assign values of E to E_ext and H to H_ext according to n
E_ext(n) = E;
H_ext(n) = H;


%% Project quantities of H on new edges in H_ext

% Project parts of Hx and Hz on boundary y = 0
n1 = 1 + ((1:nx)-1) * Mx;
n2 = 1 + ((2:nx_ext-1)-1) * Mx_ext;
H_ext(n2) = H(n1);
H_ext(n2+2*np_ext) = H(n1+2*np);

% Project parts of Hx and Hz on boundary y = y_max
n1 = 1 + ((1:nx)-1) * Mx + ((ny-1)-1) * My;
n2 = 1 + ((2:nx_ext-1)-1) * Mx_ext + ((ny_ext-1)-1) * My_ext;
H_ext(n2) = H(n1);
H_ext(n2+2*np_ext) = H(n1+2*np);

% Project parts of Hy on boundary x = 0
n1 = 1 + ((1:ny)-1) * My;
n2 = 1 + ((2:ny_ext-1)-1) * My_ext;
H_ext(n2+np_ext) = H(n1+np);
H_ext(n2+2*np_ext) = H(n1+2*np);

% Project parts of Hy on boundary x = x_max
n1 = 1 + (nx-1-1) * Mx + ((1:ny)-1) * My;
n2 = 1 + (nx_ext-1-1) * Mx + ((2:ny_ext-1)-1) * My_ext;
H_ext(n2+np_ext) = H(n1+np);
H_ext(n2+2*np_ext) = H(n1+2*np);


%% Calculate power through each relevant surface in y-direction 

% Calculate z-component of the electric field at each surface
Ez = 1/2 * (E_ext(n_rel+2*np_ext) + E_ext(n_rel+2*np_ext+Mx_ext));

% Calculate x-component of the magnetic field at each surface
Hx = 1/4 * (H_ext(n_rel) + H_ext(n_rel+Mx_ext) + H_ext(n_rel-My_ext) + H_ext(n_rel-My_ext+Mx_ext));

% Calculate x-component of the electric field at each surface
Ex = E_ext(n_rel);

% Calculate z-component of the magnetic field at each surface
Hz = 1/2 * (H_ext(n_rel+2*np_ext) + H_ext(n_rel+2*np_ext-My_ext));

% Calculate Poyntinvector to each surface
Sy = Ez .* conj(Hx) - Ex .* conj(Hz);

% Set entries at ghost edges to zero
ghost = 1 + (nx-1) * Mx + ((1:ny)-1) * My;
Sy(ghost) = 0;


%% Calculate power through each relevant surface in x-direction 

% Calculate y-component of the electric field at each surface
Ey = E_ext(n_rel+np_ext);

% Calculate z-component of the magnetic field at each surface
Hz = 1/2 * (H_ext(n_rel+2*np_ext) + H_ext(n_rel+2*np_ext-Mx_ext));

% Calculate z-component of the electric field at each surface
Ez = 1/2 * (E_ext(n_rel+2*np_ext) + E_ext(n_rel+2*np_ext+My_ext));

% Calculate y-component of the magnetic field at each surface
Hy = 1/4 * (H_ext(n_rel+np_ext) + H_ext(n_rel+np_ext-Mx_ext) + H_ext(n_rel+np_ext+My_ext) + H_ext(n_rel+np_ext+My_ext-Mx_ext));

% Calculate Poyntinvector to each surface
Sx = Ey .* conj(Hz) - Ez .* conj(Hy);

% Set entries at ghost edges to zero
ghost = 1 + ((1:nx)-1) * Mx + (ny-1) * My;
Sx(ghost) = 0;


%% Set up complete Poyntinvector and apply TD/FD

% Set up result vector
S = [Sy; Sx; zeros(np,1)];

% Multiply with 0.5, if imaginary
if isreal(Sy) == false
    S = 1/2 *S;
end

% Multiply entries in S with length of corresponding surface in the
% xy-plane
S = da * S;





