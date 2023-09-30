function [bc, W, e_exi, jsbow] = apply_bc(msh, bc, e_exitation, jsbow_excitation)
% apply_bc realizes the choosen boundary conditions, stored in bc and 
% generates the bc struct, which contains important information for the
% solvers to realize the corresponding bc as well as the projection matrix
% W, in order to realize the homogenious boundary conditions PEC or PMC
%
%          L3                 ^
%       +++++++               |
%   L4  +     +  L2           | 
%       +     +             y | 
%       +++++++               |   x
%          L1                 o------->
%
% Inputs:
%   msh                 - Mesh struct
%   bc                  - Boundary condition struct, including the boundary 
%                         definition by the user through e.g. 
%                         bc.bc = ["PMC", "OPEN", "PMC", "PEC"].
%                         The order of the conditions is: 
%                         ["L1", "L2", "L3", "L4"]
%                         Another possible input is the number of PML-cells 
%                         on each side through e.g. bc.NPML = [20,20,20,20]
%                         ONLY bc.bc is necessary, if NPML is not defined 
%                         apply std. value for number of cells -> 20.
%                         If no PML desired -> 0 entry for associated side
%   e_exitation         - Vector for electric field excitation with size 
%                         (3*np x 1) with NaNs at the index of each edge,
%                         where no excitation is defined -> Only entries
%                         unequal to NaN at indices associated with
%                         excitations with the value of the excitation in
%                         frequency domain and 1 in time domain
%   jsbow_excitation    - Vector for current excitation. Same principle as
%                         for e_exitation
%
% Outputs:
%   bc                  - Boundary condition struct with information for
%                         the solver functions regarding boundary
%                         conditions
%   W                   - Projection matrix W size (3*np x number DoFs)
%   e_exi               - New vector for electric field excitation for the
%                         solver functions size (3*np x 1)
%   jsbow               - New vector for current excitation for the
%                         solver functions size (3*np x 1)


% Get mesh parameters
Mx = msh.Mx;
My = msh.My;
nx = msh.nx;
ny = msh.ny;
np = msh.np;


%% Analyse given boundary conditions

% Check if given boundary conditions are correct
if isempty(setdiff(bc.bc, ["PMC", "PEC", "OPEN"])) == false
    error("At least one of the given boundary conditions in bc.bc is not defined or spelled wrong!")
end

% Check for open boundary conditions
open_bc = ismember(bc.bc,"OPEN"); 
bc.open_bc = open_bc; 

% In case open boundary conditions are defined, PML could be needed.
% Check if values for number of PML edges are given, else apply standard 
% value in case FD is used
if any(ismember(fieldnames(bc),"NPML"))
    bc.NPML = bc.NPML .* open_bc;
else 
    bc.NPML = [20,20,20,20] .* open_bc;
end


%% Create excitation vectors jsbow and e_exi for the solver functio

idx_exci_j = find(~isnan(jsbow_excitation));
jsbow = sparse(3*np,1);
jsbow(idx_exci_j) = jsbow_excitation(idx_exci_j);

idx_exci_e = find(~isnan(e_exitation));
e_exi = sparse(3*np,1);
e_exi(idx_exci_e) = e_exitation(idx_exci_e);


%% Create projector matrix W

% Get DoFs, depending on electric or magnetic boundary conditions and
% create the vector P with size 3*np, which has only ones at inidices for
% DoFs
P = ones(3*np,1);

% Get all edges relevant for the calculation and first set these indices as
% DoFs
idxGhostEdges = getGhostEdges_2D(msh);
P(idxGhostEdges) = 0;

% Set excitation indices to zero
P(idx_exci_e) = 0; 

% Delete indices of edges, which must be set to zero for electric boundary
% conditions

% Electric boundary for y=ymin
if strcmp(bc.bc(1), "PEC")
    n = 1 + ([1:nx]-1) * Mx; %#ok<NBRAK1> 
    P(n) = 0;
    P(n+2*np) = 0;
end

% Electric boundary for x=xmax
if strcmp(bc.bc(2), "PEC")
    n = 1 + (nx-1) * Mx + ([1:ny]-1) * My + np; %#ok<NBRAK1> 
    P(n) = 0;
    P(n+np) = 0;
end

% Electric boundary for y=ymax
if strcmp(bc.bc(3), "PEC")
    n = 1 + ([1:nx]-1) * Mx + (ny-1) * My; %#ok<NBRAK1> 
    P(n) = 0;
    P(n+2*np) = 0;
end

% Electric boundary for x=xmin
if strcmp(bc.bc(4), "PEC")
    n = 1 + ([1:ny]-1) * My + np; %#ok<NBRAK1> 
    P(n) = 0;
    P(n+np) = 0;
end

% Extract indices of DoFs from P (Indices of nonzeros)
idx_dof = find(P);

% Create P matrix -> Diagonal with ones at indices of DoFs, else 0
% Size 3*np x 3*np
P = spdiags(P, 0, 3*np, 3*np);

% Create Q matrix -> Size (3*np x number DoFs) with ones in each row in
% column i for the i-th DoF
Q = spdiags(ones(3*np,1), 0, 3*np, 3*np);
Q = Q(:,idx_dof); %#ok<FNDSB> 

% Determine projector W
W = P*Q;


end
