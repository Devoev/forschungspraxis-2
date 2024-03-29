function [ebow, hbow] = solve_helmholtz_2d_fd(msh, W, c, meps, mmui, jsbow, e_exi, f, bc)
% solve_helmholtz_2d_fd solves the 2D Helmholtz equation with a 
% current and a electric field excitation in frequency domain.
%
% Inputs:
%   msh     - Mesh struct
%   W       - Projector matrix size (3*np x number DoFs)
%   c       - Primary curl matrix
%   meps    - Permittivity matrix
%   mmui    - Inverse permeability matrix
%   jsbow   - Current excitation vector size (3*np x 1)
%             generated by function apply_bc
%   e_exi   - Electric field excitation vector (excitation voltages)
%             generated by function apply_bc;
%             Unit in x- and y-direction: V; Unit in z-direction: V/m;
%             Size (3*np, 1)
%   f       - Frequency for all quantities
%   bc      - Boundary condition struct generated by function apply_bc
%
% Outputs:
%   ebow    - Integrated electric field
%             Unit in x- and y-direction: V; Unit in z-direction: V/m;
%             Size (3*np, 1)
%   hbow    - Integrated magnetic field
%             Unit in x- and y-direction: A; Unit in z-direction: A/m;
%             Size (3*np, 1)


% Calculate omega from excitation frequency
omega = 2 * pi * f;

% add PML condition to material matrices
[meps, mmui] = calcpml_2D(msh, bc, meps, mmui);

% System matrix and rhs
A = W' * (omega^2 * meps - c' * mmui * c) * W;
rhs = W' * (1j * omega * jsbow - (omega^2 * meps - c' * mmui * c) * e_exi);

% solve system of equation equation
% tol = 1e-4;
% maxit = 600;
% e_DoF = gmres(A,rhs,[],tol,maxit);
e_DoF = A\rhs;
ebow = e_exi + W * e_DoF;

% Post processing
bbow = -c*ebow / (1i*omega);
hbow = mmui*bbow;


end
