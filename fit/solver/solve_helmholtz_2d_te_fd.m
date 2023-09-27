function [ebow, hbow] = solve_helmholtz_2d_te_fd(msh, eps, mui, jsbow, ebow_bc, omega, npml)
% SOLVE_HELMHOLTZ_2D_TE_FD Solves the 2D Helmholtz equation in the TE case in frequency domain.
%
% Inputs:
%   msh     - Mesh struct.
%   eps     - Permittivity values.
%   mui     - Reluctivity values.
%   jsbow   - Integrated current excitation.
%   ebow_bc - Boundary values of integrated electric field. NaN for DOF values.
%   omega   - Radial excitation frequency.
%   npml    - PML boundary conditions. Array of length 4.
%
% Outputs:
%   ebow    - Integrated electric field.
%   hbow    - Integrated magnetic field.

    % UPML setting
    NGRID = [msh.nx, msh.ny];

    [c, ~, ~] = createTopMatsTE(msh);

    [ds, dst, da, dat] = createGeoMatsTE(msh);

    meps = createMepsTE(msh, ds, da, dat, eps);
    mmui = createMmuiTE(msh, ds, dst, da, mui);

    % UPML tensoren - TM mode
    [sx, sy] = calcpml2D(NGRID, npml);
    sx_v = reshape(sx', [], 1);
    sy_v = reshape(sy', [], 1);
    s_mmui = sparse(diag([sy_v./sx_v; sx_v./sy_v]));
    s_eps = sparse(diag(sx_v.*sy_v));

    % UPML material matrices
    meps = meps*s_eps; % TODO: Move to bc in createMeps/createMmui
    mmui = mmui*s_mmui;

    % System matrix and rhs
    A = -c'*mmui*c + omega^2*meps;
    b = 1j*omega*jsbow;

    % Deflate system matrix
    idx_dof = isnan(ebow_bc);
    idx_bc = ~idx_dof;
    ebow_bc = ebow_bc(idx_bc);
    b = b(idx_dof) - A(idx_dof, idx_bc) * ebow_bc;
    A = A(idx_dof, idx_dof);

    % solve equation
    ebow = sparse(msh.np, 1);
    ebow(idx_dof) = A\b;
    ebow(idx_bc) = ebow_bc;

    % Post processing
    bbow = -c*ebow / (1i*omega);
    hbow = mmui*bbow;
end