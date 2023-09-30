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
xmesh = linspace(0,2,300);
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

% Edit boundary conditions
bc.bc = ["PMC", "OPEN", "PMC", "OPEN"];
bc.NPML = [30, 40, 30, 40];

%% Create Excitation

% Frequency for harmonic excitation
f = 5e9;

% Time parameters
dt = 1e-11;
tend = 5.75/f;
steps = ceil(tend/dt);

% Create empty jsbow_space vector
jsbow_excitation = NaN(3*np, 1);

% Create excitation vector for the electric field
e_exitation = NaN(3*np, 1);

% Determine indices for electric field excitation
x_L = ceil(msh.nx/2);
y_L = ceil(msh.ny/2);
n = 1 + (x_L-1)*Mx + ([1:ny]-1)*My + 2*np; %#ok<NBRAK1> 

% Set corresponding entries in e_exi to 1 
% -> Use corresponding edges as source
e_exitation(n) = 1j;  




%% Generate matrices for calculation

% Create curl, source and geometric matirces
[c, s, st] = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);

% Create permittivity matrix and it's inverse
rel_eps = boxMesher_2D(msh, boxesEps, eps0);
meps = createMeps_2D(msh, ds, da, dat, rel_eps, eps0);
mepsi = nullInv(meps);

% Create permeability matrix and it's inverse
rel_mui = boxMesher_2D(msh, boxesMui, mu0i);
mmui = createMmui_2D(msh, ds, dst, da, rel_mui, mu0i);
mmu = nullInv(mmui);



%% Apply boundary conditions and get excitation vectors for the simulation
[bc, W, e_exi, jsbow] = apply_bc(msh, bc, e_exitation, jsbow_excitation);


%% Simulate
[ebow, hbow] = solve_helmholtz_2d_fd(msh, W, c, meps, mmui, jsbow, e_exi, f, bc);


%% Plot solution for electric field in z-direction

figure(1)
idx2plot = 2*np+1:3*np;
[X,Y] = meshgrid(msh.xmesh, msh.ymesh);
e_surf = reshape(real(ebow(idx2plot)), [msh.nx, msh.ny]);
e_surf_plot = surf(X,Y,e_surf');
set(e_surf_plot,'LineStyle','none')
set(gca,'ColorScale','log')
drawnow
