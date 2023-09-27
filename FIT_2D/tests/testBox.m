%% Paths

% Clear all
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


%% Edit calculation domain

% Parameters for grid
xmesh = linspace(0,3,4);
ymesh = linspace(0,3,4);

% Parameters for permittivity
eps0 = 8.854e-12;
boxesEps(1).box = [1, 4, 1, 2];
boxesEps(1).value = 2;
boxesEps(2).box = [1, 4, 2, 4];
boxesEps(2).value = 1;

% Parameters for inverse permeability
mu0i = 1/(pi*4e-7);
boxesMui(1).box = [1, 4, 1, 2];
boxesMui(1).value = 1/2;
boxesMui(2).box = [1, 4, 2, 4];
boxesMui(2).value = 1;


%% Create mesh

% Create basic mesh object
msh = cartMesh_2D(xmesh, ymesh);

% Get geometrical matrices
[ds, dst, da, dat] = createGeoMats_2D(msh);

% Create vector with distribution of relative permittivity according to 
% boxes
rel_eps = boxMesher_2D(msh, boxesEps, eps0);

% Create permittivity matrix
meps = createMeps_2D(msh, ds, da, dat, rel_eps, eps0);

% Create vector with distribution of relative permeability according to 
% boxes
rel_mui = boxMesher_2D(msh, boxesMui, mu0i);

% Create inverse permeability matrix
mmui = createMmui_2D(msh, ds, dst, da, rel_mui, mu0i);






