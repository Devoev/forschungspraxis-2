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
xmesh = linspace(0,1,111);
ymesh = linspace(0,1,111);

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
boxesKaps(1).value = 0.03;


%% Edit excitation

% Current
I = 1;

% Frequency 
f = 2e9;


%% Generate mesh and matrices for calculation

% Create curl, source and geometric matirces
[c, s, st] = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);

% Create permittivity matrix and it's inverse
rel_eps = boxMesher_2D(msh, boxesEps, eps0);
Meps = createMeps_2D(msh, ds, da, dat, rel_eps, eps0);
Mepsi = nullInv(Meps);

% Create conductivity matrix and it's inverse
kaps_vec = boxMesher_2D(msh, boxesKaps, kap0);
Mkaps = createMeps_2D(msh, ds, da, dat, kaps_vec, 1);
Mkapsi = nullInv(Mkaps);

% Create permeability matrix and it's inverse
rel_mui = boxMesher_2D(msh, boxesMui, mu0i);
Mmui = createMmui_2D(msh, ds, dst, da, rel_mui, mu0i);
Mmu = nullInv(Mmui);


%% Create excitation vector

% Create empty jsbow_space vector
jsbow = sparse(3*msh.np, 1);

% Determine index in the middle of the calculation domain
x_L = ceil(msh.nx/2);
y_L = ceil(msh.ny/2);
n = 1 + x_L*Mx + y_L*My + 2*np;

% Set current on the determined index
jsbow(n) = I;  
omega = 2 * pi * f;

% Set up stystem of equations
A = -c'*Mmui*c + omega^2*Meps - 1i*omega*Mkaps;
b = 1j*omega*jsbow;

% Deflate system matrix
idx_dof = getNotGhostEdges_2D(msh);
idx_dof = nonzeros(idx_dof);
b = b(idx_dof);
A = A(idx_dof, idx_dof);

% Solve system
ebow = sparse(3*msh.np, 1);
ebow(idx_dof) = A\b;

% Plot solution
ebow = real(ebow);

zlimit = 700;

z_plane = 1;
idx2plot = 2*np+1:3*np;
ebow_mat = reshape(ebow(idx2plot),nx,ny);
figure(1)
mesh(ebow_mat)
xlabel('i')
ylabel('j')
zlabel(['z-Komponente des E-Feldes für z=',num2str(z_plane)])
axis([1 nx 1 ny -zlimit zlimit])

drawnow