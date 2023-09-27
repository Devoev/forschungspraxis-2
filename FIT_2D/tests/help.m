
%% Edit calculation domain
clc
clear

% Steps of mesh
xmesh = linspace(0,1,500);
ymesh = linspace(0,1,500);

% Create basic mesh object
msh = cartMesh_2D(xmesh, ymesh); 
Mx = msh.Mx;
My = msh.My;
nx = msh.nx;
ny = msh.ny;
np = msh.np;

% Parameters for permittivity
eps0 = 8.854e-12;
boxesEps(1).box = [1, nx, 1, ny];
boxesEps(1).value = 1;

% Parameters for inverse permeability
mu0i = 1/(pi*4e-7);
boxesMui(1).box = [1, nx, 1, ny];
boxesMui(1).value = 1;

% Relative permeability and permittivity
mui = 1;
epsilon = 1;


%% Add paths

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


%% Generate mesh and matrices for calculation

% Create curl, source and geometric matirces
[c, s, st] = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);

% Create permittivity matrix and it's inverse
rel_eps = boxMesher_2D(msh, boxesEps, eps0);
meps = createMeps_2D(msh, ds, da, dat, rel_eps, eps0);
Mepsi = nullInv(meps);

% Create permeability matrix and it's inverse
rel_mui = boxMesher_2D(msh, boxesMui, mu0i);
mmui = createMmui_2D(msh, ds, dst, da, rel_mui, mu0i);
Mmu = nullInv(mmui);


%% PML

npml = [20, 20, 20, 20];

NGRID = [msh.nx, msh.ny];

% UPML tensoren - TM mode
[sx, sy] = calcpml_2D(NGRID, npml);
sx_v = reshape(sx', [], 1);
sy_v = reshape(sy', [], 1);
s_mmui = sparse(diag([sy_v./sx_v; sx_v./sy_v]));
s_eps = sparse(diag(sx_v.*sy_v));

% UPML material matrices
meps = repmat(s_eps,3,3); 
mmui = mmui*s_mmui;

s_eps2 = spdiags(s_eps,3*np,3*np);


