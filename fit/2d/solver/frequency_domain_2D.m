function [ebow, hbow] = frequency_domain_2D(msh, c, meps, mmui, mkaps, jsbow, idx_bc, ebow_bc, omega, NPML)
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
%   bc      - Boundary conditions.
%
% Outputs:
%   ebow    - Integrated electric field.
%   hbow    - Integrated magnetic field.

    % add PML condition to material mats
    [meps, mmui] = calcpml_2D(msh, NPML, meps, mmui);

    % System matrix and rhs

    % If there is a empty kappa matrix
    if isempty(nonzeros(mkaps)) 
 
        A = -c'*mmui*c + omega^2*meps;
        b = 1j*omega*jsbow;

    % If a conductivity is given
    else

        A = -c'*mmui*c + omega^2*meps - 1i*omega*mkaps;
        b = 1j*omega*jsbow;
        
    end 

    % Deflate system matrix
    idx_dof = nonzeros(getNotGhostEdges_2D(msh));
    b = b - A(:, idx_bc) * ebow_bc;
    idx_dof = setdiff(idx_dof, idx_bc);
    b = b(idx_dof);
    A = A(idx_dof, idx_dof);

    % solve equation
    ebow = sparse(msh.np, 1);
    ebow(idx_dof) = A\b;
    ebow(idx_bc) = ebow_bc;

    % Post processing
    bbow = -c*ebow / (1i*omega);
    hbow = mmui*bbow;

end
