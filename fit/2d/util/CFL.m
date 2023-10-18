function tmax = CFL(msh, MAT)
% CFL calculates the maximum timestep applyable to the Leapfrog algorithmn
% according to the CFL criterion
%
% Input
% msh               - mesh struct
% MAT               - MAT struct containing matrices for calculations in
%                     FIT
%
% Output
% tmax              - maximum time step according to CFL criterion

% Get minimal edge lengths in all direction
dx = min(diff(msh.xmesh));
dy = min(diff(msh.ymesh));
dz = msh.lz;

% Calculate maximum time step according to CFL
tmax = sqrt(MAT.epsilon0/MAT.mu0i) * sqrt(1/(1/(dx^2) + 1/(dy^2) + 1/(dz^2)));

end
