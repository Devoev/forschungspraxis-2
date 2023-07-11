function [sx, sy] = calcpml2D(NGRID, NPML)
% calcpml2d. von A.G.
% 
% Inputs:
%   NGRID   - mesh dimentions list with [x-dim, y-dim]
%   NPML    - number of mesh layers in each direction [x-layers, y-layers]
% 
% Outputs:
%   sx    - x pml tensor w. dimension like mesh Nx x Ny
%   sy    - y pml tensor w. dimension like mesh Nx x Ny

% pml settings
sigma_max = 1;
a_max = 3;
p = 3;

% input
Nx = NGRID(1);
Ny = NGRID(2);
NXHI = NPML(1); % +X
NYLO = NPML(2); % +Y
NXLO = NXHI; % -X
NYHI = NYLO; % -Y
Lx = NXHI;
Ly = NYLO;

% init tensors ------------------------------------------------------------
imp0 = 376.73;
sx = sparse(ones(Nx, Ny));
sy = sparse(ones(Nx, Ny));

% computing pml params
sigmax =@(x) sigma_max*sin(pi*x/(2*Lx))^2;
sigmay =@(y) sigma_max*sin(pi*y/(2*Ly))^2;

ax =@(x) 1+a_max*(x/Lx)^p;
ay =@(y) 1+a_max*(y/Ly)^p;

% building param functions
f_sx =@(x) ax(x)*(1+1j*imp0*sigmax(x));
f_sy =@(y) ay(y)*(1+1j*imp0*sigmay(y));

% calc sx tensor ----------------------------------------------------------
% add +X direction PML
for nx = 1:NXHI
    val = nx/NXHI;
    % set column from inside to outside
    sx(Nx-NXHI+nx,:) = f_sx(val); 
end
% add -X direction PML
for nx = 1:NXLO
    val = nx/NXLO;
    % set column from inside to outside
    sx(NXLO-nx+1,:) = f_sx(val); 
end

% calc sy tensor ----------------------------------------------------------
% add +Y direction PML
for ny = 1:NYLO
    val = ny/NYLO;
    % set row from inside to outside of boundary
    sy(:,NYLO-ny+1) = f_sy(val);
end
% add -Y direction PML
for ny = 1:NYHI
    val = ny/NXHI;
    % set row from inside to outside of boundary
    sy(:,Ny-NYHI+ny) = f_sy(val);
end

end