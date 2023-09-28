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

% Set open boundary for each side (false -> PMC)
open_bc = [true, true, true, true];  % [L1, L2, L3, L4];

% Relative permeability and permittivity
mui = 1;
epsilon = 1;


%% Edit excitation

% Gauss
sigma = 6e-10;

% Harmonic 
f = 1e9;

% Time parameters
dt = 1e-11;
tend = 10*sigma;
steps = floor(tend/dt);

% Edit sourcetype
sourcetype= 2;  % 1: Gauss Anregung, 2: Harmonisch, 3: Konstante Anregung


%% Generate mesh and matrices for calculation

% Create curl, source and geometric matirces
[c, s, st] = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);

% Create permittivity matrix and it's inverse
rel_eps = boxMesher_2D(msh, boxesEps, eps0);
Meps = createMeps_2D(msh, ds, da, dat, rel_eps, eps0);
Mepsi = nullInv(Meps);

% Create permeability matrix and it's inverse
rel_mui = boxMesher_2D(msh, boxesMui, mu0i);
Mmui = createMmui_2D(msh, ds, dst, da, rel_mui, mu0i);
Mmu = nullInv(Mmui);

% Determine z-edges for open boundary
[mur_edges, mur_n_edges, mur_deltas] = initMur_2D(msh, open_bc);


%% Create excitation vector

% Create empty jsbow_space vector
jsbow_space = zeros(3*np, 1);

% Determine index in the middle of the calculation domain
x_L = ceil(msh.nx/2);
y_L = ceil(msh.ny/2);
n = 1 + x_L*Mx + y_L*My + 2*np;

% Set current on the determined index
jsbow_space(n) = 1;  
jmax = 1;

% Gauss excitation
jsbow_gauss = @(t)(jsbow_space * jmax * exp(-4*((t-sigma)/sigma)^2));

% harmonic excitation
jsbow_harm = @(t)(jsbow_space * jmax * sin(2*pi*f*t));

% constant excitation
jsbow_const = @(t)(jsbow_space * jmax);


%% Initialize simulation

% Initialize for leapfrog
ebow_new = sparse(3*np,1);
hbow_new = sparse(3*np,1);
energy = zeros(1,steps);
power_source = zeros(1,steps);

% Plot parameter for "movie"
figure(1)
zlimit = 700;
draw_only_every = 5;


%% Execute simulation

for ii = 1:steps

    % Calculate time t
    t = ii*dt;

    % Old values
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    % Calculate value for excitation

    % Gauss
    if sourcetype == 1

        if t <= 2*sigma
            js = jsbow_gauss(t);
        else
            js = sparse(3*np,1);
        end

    % Harmonic
    elseif sourcetype == 2
        js = jsbow_harm(t);

    % Constant
    elseif sourcetype == 3
        js = jsbow_const(t);
    end
    
    % Execute timestep with leapfrog
    [hbow_new,ebow_new] = leapfrog_2D(hbow_old, ebow_old, js, Mmui, Mepsi, c, dt);
    % Apply open boundary with mur cond
    ebow_new = applyMur_2D(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt);
    
    % Draw electric field
    if mod(ii, draw_only_every)
        z_plane = 1;
        idx2plot = 2*np+1:3*np;
		ebow_mat = reshape(ebow_new(idx2plot),nx,ny);
    	figure(1)
        mesh(ebow_mat)
        xlabel('i')
        ylabel('j')
        zlabel('z-component of electric field')
		axis([1 nx 1 ny -zlimit zlimit])
		clim([-zlimit zlimit])
        drawnow
    end

    % Calculate energy in the system and power of the source
    energy_t = 0.5 * (ebow_new' * Meps*ebow_new + hbow_new' * Mmu * hbow_new);
    leistungQuelle_t = ebow_new' * js;

    % Save values
    energy(ii) =  energy_t;
    power_source(ii) = leistungQuelle_t;
end


%% Plot further results

% Plot excitation
figure(2)
jsbow_plot = zeros(1,steps);
for step = 1:steps
    if sourcetype == 1
        jsbow_spatial = jsbow_gauss(step*dt);
    elseif sourcetype == 2
        jsbow_spatial = jsbow_harm(step*dt);
    elseif sourcetype == 3
        jsbow_spatial = jsbow_const(step*dt);
    end
    nonzero_idx = find(jsbow_spatial~=0);
    jsbow_plot(step) = jsbow_spatial(nonzero_idx);
end
plot(dt:dt:dt*steps, jsbow_plot);
xlabel('t in s');
ylabel('Excitation current J in A');


% Plot energy over time
figure(3); clf;
plot (dt:dt:dt*steps, energy)
legend(['Zeitschritt: ', num2str(dt)])
xlabel('t in s')
ylabel('Energie des EM-Feldes W in J')


% Plot power
leistungSystem = diff(energy) / dt;
figure(4); clf;
hold on
plot(2*dt:dt:dt*(steps), leistungSystem)
plot(dt:dt:dt*steps, power_source, 'r')
hold off
legend('Leistung System', 'Leistung Quelle')
xlabel('t in s')
ylabel('Leistung P in W')
