function [ebow, hbow, relRes] = solveHelmholtz2D(msh, eps, mui, jsbow, omega, bc)
% SOLVE_HELMHOLTZ_2D Solves the 2D Helmholtz equation.
% Inputs:
%   msh     - Mesh struct.
%   eps     - Permittivity values.
%   mui     - Reluctivity values.
%   jsbow   - Integrated current excitation.
%   omega   - Radial excitation frequency.
%   bc      - Boundary conditions.
% Outputs:
%   ebow    - Integrated electric field.
%   hbow    - Integrated magnetic field.

    % Anzahl der Rechenpunkte des Gitters
    np = msh.np;

    [c, g, st] = createTopMats2DTE(msh);

    [ds, dst, da, dat] = createGeoMats2D(msh);

    meps = create2DMeps(msh, ds, da, dat, eps);
    msig = create2DMeps(msh, ds, da, dat, 1e-4); % TODO: Remove later
    mmui = create2DMmui(msh, ds, dst, da, mui);

    % Berechnung Systemmatrix A und rechte Seite rhs
%    A = -c'*mmui*c + omega^2*meps;
    A = -c'*mmui*c + omega^2*meps - 1j*omega*msig;
    rhs = 1j*omega*jsbow;

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