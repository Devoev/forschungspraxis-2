%% Add paths

% Clear variables
clc
clear

% Get parent directory
filePath = matlab.desktop.editor.getActiveFilename;
[ParentFolderPath] = fileparts(filePath);
parent = fileparts(ParentFolderPath) ;

% Paths to add
path_msh_func = append(parent, '\mesh');
path_mat_func = append(parent, '\matrices');
path_solver_func = append(parent, '\solver');
path_util_func = append(parent, '\util');

% Add paths
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func)


%% Edit basic calculation domain

% Steps of mesh
xmesh = linspace(0,2,300);
ymesh = linspace(0,1,200);

% Create basic mesh object
msh = cartMesh_2D(xmesh, ymesh); 
Mx = msh.Mx;
My = msh.My;
nx = msh.nx;
ny = msh.ny;
np = msh.np;

% Edit boundary conditions
bc.bc = ["PMC", "OPEN", "PMC", "OPEN"];


%% Edit material regions and add them to the object material_regions

% Add basic constants to material_regions object
material_regions.epsilon0 = 8.854187e-12;
material_regions.mu0i = 1/(pi*4e-7);

% Regions for relative permittivity
% Relative permittivity everywhere equal to one
boxesEpsilonR(1).box = [1, nx, 1, ny];
boxesEpsilonR(1).value = 1;
material_regions.boxesEpsilonR = boxesEpsilonR;

% Regions for inverse relative permeability
% Inverse relative permeability everywhere equal to one
boxesMuiR(1).box = [1, nx, 1, ny];
boxesMuiR(1).value = 1;
material_regions.boxesMuiR = boxesMuiR;


%% Create Excitation

% Frequency for harmonic excitation
f = 1e9;

% Create empty jsbow_space vector
jsbow_excitation = NaN(3*np, 1);

% Create excitation vector for the electric field
e_exitation = NaN(3*np, 1);

% Determine indices for electric field excitation
x_L = ceil(msh.nx/2);
n = 1 + (x_L-1)*Mx + ((1:ny)-1)*My + 2*np; 

% Set corresponding entries in e_exi to the value of the desired amplitude
% -> Use corresponding edges as source
e_exitation(n) = 1;  


%% Apply boundary conditions and get excitation vectors for the simulation

[bc, W, e_exi, jsbow] = apply_bc(msh, bc, e_exitation, jsbow_excitation);


%% Generate matrices for calculation

[MAT, bc] = generate_MAT(msh, bc, material_regions, ["CurlP"]); %#ok<NBRAK2> 


%% Simulate

[ebow, hbow] = solve_helmholtz_2d_fd(msh, W, MAT.c, MAT.meps, MAT.mmui, jsbow, e_exi, f, bc);


%% Plot solution for electric field in z-direction

figure(1)
idx2plot = 2*np+1:3*np;
[X,Y] = meshgrid(msh.xmesh, msh.ymesh);
e_surf = reshape(real(ebow(idx2plot)), [msh.nx, msh.ny]);
e_surf_plot = surf(X,Y,e_surf');
set(e_surf_plot,'LineStyle','none')
set(gca,'ColorScale','log')
drawnow
