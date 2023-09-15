function [ebow, hbow, relRes] = solveHelmholtzTE(msh, eps, mui, jsbow, ebow_bc, omega, npml)
% SOLVE_HELMHOLTZ_TE Solves the 2D Helmholtz equation in the TE case.
%
% Inputs:
%   msh     - Mesh struct.
%   eps     - Permittivity values.
%   mui     - Reluctivity values.
%   jsbow   - Integrated current excitation.
%   ebow_bc - Boundary values of integrated electric field.
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
    rhs = 1j*omega*jsbow;

    % Deflate/ Inflate system matrix
    % TODO: Inflation with ebow_bc
    idx_dof = isnan(ebow_bc);
    idx_bc = ~idx_dof;

    % solve equation
%    [ebow, flag, relRes, iter, resVec] = gmres(A, rhs, 20, 1e-10, 1000); % TODO: direct vs iteratve?
%    if flag == 0
%      fprintf('gmres(20): converged at iteration %2d to a solution with relative residual %d.\n',iter,relRes);
%    else
%      error('gmres(20): some error ocurred, please check flag output.')
%    end
%    relRes = resVec./norm(rhs);
    ebow = A\rhs;

    % Post processing
    bbow = -c*ebow / (1i*omega);
    hbow = mmui*bbow;
end