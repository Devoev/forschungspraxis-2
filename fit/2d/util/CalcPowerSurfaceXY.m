function [S] = CalcPowerSurfaceXY(msh, ebow, hbow, ds, dst, da)
% CalcPowerSurfaceXY calculates the power through each surface in x- and
% y-direction utilizing the Poyntinvector. Applyable for time and frequency
% domain. !!! The power is only calculated on surfaces on a subdomain,
% which is one point shorter in each direction. Only faces inside or on the
% boundary of the subdomain are considered !!!
%
% Inputs:
%   msh                 - Mesh struct
%   ebow                - Vector with integrated electric field
%   hbow                - Vector with integrated magnetic field
%   ds                  - Primary edge matrix
%   dst                 - Dual edge matrix
%   da                  - Primary surface matrix
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


%% Calculate indices of all points in the subdomain relevant for surfaces 
% in y-direction
n_subdom = 1 + ((2:nx-2)-1) * Mx;
n_subdom = repmat(n_subdom, ny-2, 1)' + My * (1:ny-2);
n_subdom = reshape(n_subdom, (nx-3)*(ny-2), 1);


%% Calculate power through each relevant surface in y-direction 

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


%% Calculate indices of all points in the subdomain relevant for surfaces 
% in x-direction
n_subdom = 1 + ((2:nx-1)-1) * Mx;
n_subdom = repmat(n_subdom, ny-3, 1)' + My * (1:ny-3);
n_subdom = reshape(n_subdom, (nx-2)*(ny-3), 1);


%% Calculate power through each relevant surface in x-direction 

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
