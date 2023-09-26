%% Task3: Thin-film interference

clc
clear all

% Path mesh functions
path_msh_func = './fit/mesh';
path_mat_func = './fit/matrices';
path_solver_func = './fit/solver';
path_util_func = './fit/util';
path_verify_func = './task1/verifications';

% Add paths
cd('../');
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func, path_verify_func)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_field = 1;


%% Problem Definition
c = 3e8;            % [m/s]
eps = 8.854e-12;
mui = 1/(4*pi*1e-7);
n = 2;

v0 = c;             % [m/s]
v1 = c/n;           % [m/s] given c = 1/sqrt(mu0 * eps0)
% refraction idx n=2 by multiplying mu0 or eps0 for thin film area by 4
eps1 = 4*eps;
% TODO: inhomogenous Meps matrix

lambda1 = 430e-9;   % [m]
f1 = c/lambda1;     % [Hz]
omega1 = 2*pi*f1;   % [1/s]
E1 = 250;           % [V/m]

lambda2 = 510e-9;   % [m]
f2 = c/lambda2;     % [Hz]
omega2 = 2*pi*f2;   % [1/s]
E2 = 500;           % [V/m]

%         thin film geometry           |   
%        -h/2    L2[PML]  h/2          |      b: middle index
%         |------ y ------|            |         
%        # | | | | | | | | #           |      d: distance between the
%  [PEC] # v v v v v v v v #   [PEC]   |         excitations in wave-
%  L3    #                 #   L1      |         length  
%  x     #-----------------#border1(L/2)|
%        #----thin film----#  a        |      L: side index bc
% border2#-----------------#           |         [L1, L2, L3, L4]
%        #                 #           |
%        ###################           |
%               L4 [PML]

%% geometry
h = 4e-6;   % [m]
L = 10e-6;  % [m]
a = 100e-9; % [m]

NPML = [1, 20, 1, 20];  % [L1, L2, L3, L4]; 0,1:=PMC

%% mesh
elem_per_wavelength = 15;
dx = lambda1*(NPML(2)+NPML(4))/elem_per_wavelength;  % Extra space in +x direction
xmesh = linspace(0, L + 2*dx, ceil( (L + 2*dx)/lambda1*elem_per_wavelength) );
ymesh = linspace(-h/2, h/2, ceil( h/lambda1*elem_per_wavelength ));
msh = cartMesh2D(xmesh, ymesh);

% border x indices of different permittivity == thin film borders
border1_x_idx = round((L/2)/(L/msh.nx)); 
border2_x_idx = border1_x_idx + round(a/(L/msh.nx));
actual_thickness = (border2_x_idx-border1_x_idx)*(L/msh.nx)

%% excitation
indices = 1:msh.np;
idx = indices(~(mod(indices,msh.nx)-31));  % idx of all L2 boundary elements

jsbow = sparse(msh.np, 1);
ebow1_bc = NaN(msh.np, 1);
ebow2_bc = NaN(msh.np, 1);
ebow1_bc(idx) = E1;
ebow2_bc(idx) = E2;
% Anregung erstmal in z-Richtung (wie Task1). entspricht Aufgabe d)

%% solve system
ebow1 = solveHelmholtzTE(msh, eps, mui, jsbow, ebow1_bc, omega1, NPML);
ebow2 = solveHelmholtzTE(msh, eps, mui, jsbow, ebow2_bc, omega2, NPML);
ebow = ebow1 + ebow2;

%% Postprocessing
if plot_field
    figure
    [X,Y] = meshgrid(msh.xmesh, msh.ymesh);
    e_surf = reshape(real(ebow), [msh.nx, msh.ny]);
    e_surf_plot = surf(X,Y,e_surf');
    xlabel('X');
    ylabel('Y');
    zlim([-E2*6 E2*6])
    set(e_surf_plot,'LineStyle','none')
    set(gca,'ColorScale','log')
end


