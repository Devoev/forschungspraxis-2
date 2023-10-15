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

polarization = 'y';
excitation = 3; % 1: x-transversal wave, 2: y-transversal wave
                % 3: point in middle of domain

% Steps of mesh
xmesh = linspace(0,2,200);
ymesh = linspace(0,1,100);

% Create basic mesh object
msh = cartMesh_2D(xmesh, ymesh); 
Mx = msh.Mx;
My = msh.My;
nx = msh.nx;
ny = msh.ny;
np = msh.np;

% Edit boundary conditions
bc.bc = ["OPEN", "OPEN", "OPEN", "OPEN"];


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

% Harmonic 
f = 5e8;

% Time parameters
dt = 1e-11;
tend = 4.5/f;
steps = ceil(tend/dt);

% Create empty jsbow_space vector
jsbow_excitation = NaN(3*np, 1);

% Create excitation vector for the electric field
e_exitation = NaN(3*np, 1);

% Determine indices for electric field excitation
x_h = ceil(msh.nx*0.5);
y_h = ceil(msh.ny*0.5);

if polarization == 'x'
    pol_offset = 0;
    idx2plot = 1:msh.np;
elseif polarization == 'y'
    pol_offset = msh.np;
    idx2plot = msh.np+1:2*msh.np;
else
    pol_offset = 2*msh.np;
    idx2plot = 2*msh.np+1:3*msh.np;
end

if excitation == 1
    n = 1 + (x_h-1)*Mx + ((1:ny)-1)*My + pol_offset;
elseif excitation == 2
    n = 1 + (y_h-1)*My + (1:nx)-1 + pol_offset;
else
    n = 1 + (x_h-1)*Mx + (y_h-1)*My +pol_offset;
end
% Set corresponding entries in e_exi to one, so these edges can be used as
% sources
e_exitation(n) = 1;  


%% Apply boundary conditions and get excitation vectors for the simulation

[bc, W, e_exi, jsbow] = apply_bc(msh, bc, e_exitation, jsbow_excitation);


%% Generate matrices for calculation

[MAT, bc] = generate_MAT(msh, bc, material_regions, ["CurlP"]); %#ok<NBRAK2> 


%% Initialize simulation in time domain with leapfrog

% Function for excitation
s_source = e_exi;
e_harm = @(t)(s_source * 1 * sin(2*pi*f*t));

% Initialize open boundary condition if needed
[mur_edges,mur_n_edges, mur_deltas] = initMur_2D(msh, bc);

% Initialize ebow and hbow
ebow_new = zeros(3*np,1);
hbow_new = zeros(3*np,1);

% Add inverse permittivity matrix
MAT.mepsi = nullInv(MAT.meps);

% Plot parameter for "movie"
figure(1)
zlimit = 2.5;
draw_only_every = 10;


%% Execute simulation by applying Leapfrog algorithm
for ii = 0:steps-1

    % Calculate time t
    t = ii*dt;

    % Old values
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    % Calculate value for excitation
    e_exi_old = e_harm(t);
    e_exi_new = e_harm(t+1);

    % Execute timestep with leapfrog
    [ebow_new,hbow_new] = solve_leapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.c,dt,W);

    % Apply open boundary with mur cond
    ebow_new = applyMur_2D(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt, bc);

    % Draw electric field
    [X,Y] = meshgrid(msh.xmesh, msh.ymesh);
    if mod(ii, draw_only_every)
        f = figure(1);
        e_surf = reshape(ebow_new(idx2plot), [msh.nx, msh.ny]);
        e_surf_plot = surf(X,Y,e_surf');
        zlim([-zlimit zlimit])
        set(e_surf_plot,'LineStyle','none')
        set(gca,'ColorScale','log')
        drawnow

    end

end
