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


%% Edit calculation domain
clc
clear

% Steps of mesh
xmesh = linspace(0,1,200);
ymesh = linspace(0,1,200);

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

% Parameters for conductivity
kap0 = 0;
boxesKaps(1).box = [1, nx, 1, ny];
boxesKaps(1).value = 0;

% Distances for PML
NPML = [20,20,20,20];


%% Edit excitation

% Get z-edge in the middle of the calc domain
x_L = ceil(msh.nx/2);
y_L = ceil(msh.ny/2);
n = 1 + (x_L-1)*Mx + (y_L-1)*My + 2*np;

% Set indices for electric field excitation
idx_bc = [n]; %#ok<NBRAK2> 

% Set corresponding values for electric field excitation
ebow_bc = [250];  %#ok<NBRAK2> 

% Frequency for excitation
f = 2e9;

% Source current
jsbow = sparse(3*np,1);


%% Generate mesh and matrices for calculation

% Create curl, source and geometric matirces
[c, s, st] = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);

% Create permittivity matrix and it's inverse
rel_eps = boxMesher_2D(msh, boxesEps, eps0);
meps = createMeps_2D(msh, ds, da, dat, rel_eps, eps0);
mepsi = nullInv(meps);

% Create conductivity matrix and it's inverse
kaps_vec = boxMesher_2D(msh, boxesKaps, kap0);
mkaps = createMeps_2D(msh, ds, da, dat, kaps_vec, 1);
mkapsi = nullInv(mkaps);

% Create permeability matrix and it's inverse
rel_mui = boxMesher_2D(msh, boxesMui, mu0i);
mmui = createMmui_2D(msh, ds, dst, da, rel_mui, mu0i);
mmu = nullInv(mmui);


%% Solve in frequency domain
omega = 2 * pi * f;
[ebow, hbow] = frequency_domain_2D(msh, c, meps, mmui, mkaps, jsbow, idx_bc, ebow_bc, omega, NPML);


%% Plot solution for electric field
figure
idx2plot = 2*np+1:3*np;
[X,Y] = meshgrid(msh.xmesh, msh.ymesh);
e_surf = reshape(real(ebow(idx2plot)), [msh.nx, msh.ny]);
e_surf_plot = surf(X,Y,e_surf');
set(e_surf_plot,'LineStyle','none')
set(gca,'ColorScale','log')
drawnow
