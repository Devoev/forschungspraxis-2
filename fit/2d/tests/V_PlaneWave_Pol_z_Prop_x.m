%% Verification script

% Polarization: Z           Propagation: X

% Verification script for the simulation of a plane wave with polarization
% of the electric field in z-direction and propagation in x-direction. The
% domain is devided in two region. For x<0, the relative permittivity is 1,
% for x>0 it is 4. The simulation is executet in frequency and time domain 
% and the results are compared with an analytic solution.

% The domain is a 4 micro meter times 4 micro meter square (2D) and the
% plane wave for the exitation is placed at x = -2 micro meter.
% The wave has a wavelength of 430nm and an amplitude of 500V. 


%% Add paths

% Clear variables
clc
clear

% Get parent directory
filePath = matlab.desktop.editor.getActiveFilename;
[parent] = fileparts(filePath);
parent = fileparts(parent) ;

% Paths to add
path_msh_func = append(parent, '\mesh');
path_mat_func = append(parent, '\matrices');
path_solver_func = append(parent, '\solver');
path_util_func = append(parent, '\util');

% Add paths
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func)


%% Define excitation

% Add basic constants to material_regions object
material_regions.epsilon0 = 8.854187e-12;
material_regions.mu0i = 1/(pi*4e-7);

% Excitation 1
lambda      = 430e-9;
f           = sqrt(material_regions.mu0i/material_regions.epsilon0)/lambda;           
E           = 500;
func_exi    = @(t)(E * sin(2*pi*f*t));


%% Define important parameters for the simulation

% Elements per wavelength
elem_per_wavelength = 30;

% Offset in each direction in elements
offset = [0,25,0,25];

% Edit boundary conditions
bc.bc = ["PMC", "OPEN", "PMC", "OPEN"];


%% Edit basic calculation domain 

% Distance in x-direction in micro meters
Lx = 4;

% Distance in y-direction in micro meters
Ly = 4;

% Define maximal edge lenth in the mesh
le = lambda/elem_per_wavelength;

% Fit edge length to be an integer divisor of 0.5 micro meter and calculate
% number of edges per micro meter
num_e   = 2 * ceil(0.5e-6 / le);
le      = 1e-6 / num_e;

% Calculate number of points in x- and y-direction
points_x = num_e * Lx + 1;
points_y = num_e * Ly + 1;

% Calculate xmesh with respect to the choosen offset
x_offset1 = -Lx/2 * 1e-6  + (-offset(4):-1) * le;
x_offset2 = Lx/2 * 1e-6 + (1:offset(2)) * le;
x_basic = linspace(-Lx/2 * 1e-6, Lx/2 * 1e-6, points_x);
xmesh = [x_offset1, x_basic, x_offset2];

% Calculate ymesh with respect to the choosen offset
y_offset1 = -Ly/2 * 1e-6  + (-offset(1):-1) * le;
y_offset2 = Ly/2 * 1e-6 + (1:offset(3)) * le;
y_basic = linspace(-Ly/2 * 1e-6, Ly/2 * 1e-6, points_y);
ymesh = [y_offset1, y_basic, y_offset2];

% Create basic mesh object
msh = cartMesh_2D(xmesh, ymesh); 
Mx = msh.Mx;
My = msh.My;
nx = msh.nx;
ny = msh.ny;
np = msh.np;
lz = msh.lz;


%% Find important indices in xmesh and ymesh

% Index of x = 0
idx_x0      = find(xmesh == 0);

% Index of x = -Lx/2
idx_xmLxh   = find(xmesh == -Lx/2 * 1e-6);

% Index of x = Lx/2
idx_xpLxh   = find(xmesh == Lx/2 * 1e-6);

% Index of y = 0
idx_y0      = find(ymesh == 0);

% Index of y = -Ly/2
idx_ymLyh   = find(ymesh == -Ly/2 * 1e-6);

% Index of y = Ly/2
idx_ypLyh   = find(ymesh == Ly/2 * 1e-6);


%% Create the vectors describing the excitation

% Create empty jsbow_space vector
jsbow_excitation = NaN(3*np, 1);

% Create excitation vector for the electric field
e_exitation = NaN(3*np, 1);

% Determine indices for points at x = -Lx/2
n_exi = 1 + (idx_xmLxh - 1) * Mx + ((idx_ymLyh:idx_ypLyh) - 1) * My;

% Use corresponding edges in z-direction for excitation
e_exitation(n_exi + 2*np) = 1;


%% Edit material regions and add them to the object material_regions

% Regions for relative permittivity
% Negative x: Epsr = 1; Positive x: Epsr = 4
boxesEpsilonR(1).box = [1, idx_x0, 1, ny];
boxesEpsilonR(1).value = 1;
boxesEpsilonR(2).box = [idx_x0, nx, 1, ny];
boxesEpsilonR(2).value = 4;
material_regions.boxesEpsilonR = boxesEpsilonR;

% Regions for inverse relative permeability
% Inverse relative permeability everywhere equal to one
boxesMuiR(1).box = [1, nx, 1, ny];
boxesMuiR(1).value = 1;
material_regions.boxesMuiR = boxesMuiR;


%% Apply boundary conditions and get excitation vectors for the simulation

[bc, W, e_exi, jsbow] = apply_bc(msh, bc, e_exitation, jsbow_excitation);


%% Generate matrices for calculation

[MAT, bc] = generate_MAT(msh, bc, material_regions, ["CurlP"]); %#ok<NBRAK2> 

rel_eps = ones(np,1);
n = repmat(idx_x0:nx, ny, 1)' + (0:ny-1) * My;
n = reshape(n, 1, max(size(n) * min(size(n))));
rel_eps(n) = 4;

% Set values in ghost volumes to zero
gvx = 1 + (nx-1)*Mx + (0:ny-1) * My;
gvy = 1 + (0:nx-1)*Mx + (ny-1) * My;
rel_eps(gvx) = 0;
rel_eps(gvy) = 0;

meps = createMeps_2D(msh, MAT.ds, MAT.da, MAT.dat, rel_eps, material_regions.epsilon0);
MAT.meps = meps;


%% Set up parameters for the simulation in time domain for excitation 1 

% Time step size
dt = 0.5 * CFL(msh, MAT);

% Fit dt to period of excitation 1
dt = 1/f / ceil(1/f / dt);

% End time: 15.25 periods
t_end = floor(2e-6 * (1+2) / sqrt(material_regions.mu0i/material_regions.epsilon0) * f) /f + 1/f;


%% Simulate in time domain

% Initialize open boundary condition if needed
[mur_edges,mur_n_edges, mur_deltas] = initMur_2D(msh, bc);

% Initialize ebow and hbow
ebow_new = zeros(3*np,1);
hbow_new = zeros(3*np,1);

% Add inverse permittivity matrix
MAT.mepsi = nullInv(MAT.meps);

% Initialize calculation of avarage power
S_time = zeros(3*np,1);
i_steps = 0;

% Calculate time steps
steps_movie = 0;
for t = linspace(0,t_end,ceil(t_end/dt))

    % Old values
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    % Calculate value for excitation
    e_exi_old = e_exi * func_exi(t) * lz;
    e_exi_new = e_exi * func_exi(t+1*dt) * lz;

    % Execute timestep with leapfrog
    [ebow_new,hbow_new] = solve_leapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.c,dt,W);

    % Apply open boundary with mur condition
    ebow_new = applyMur_2D(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt, bc);

    % Sum up power through each surface for one period before the end
    if t >= t_end - 1/f
        S_time = S_time + CalcPoyntinvectorXY(msh, ebow_new, hbow_new, MAT.ds, MAT.dst);
        i_steps = i_steps + 1;
    end

    % Plot excitation of the electric field
    steps_movie = steps_movie + 1;

    if steps_movie == 10

        steps_movie = 0;

        [X,Y] = meshgrid(xmesh, ymesh);
        e_surf = reshape(ebow_new(2*np+1:3*np,1)/lz, [msh.nx, msh.ny]);
        
        figure(1)
        e_surf_plot = surf(X,Y,e_surf');
    
        xlim([-Lx/2*1e-6, Lx/2*1e-6])
        ylim([-Ly/2*1e-6, Ly/2*1e-6])
        set(e_surf_plot,'LineStyle','none')
        view(2)
        title('Amplitude of the electric field','Interpreter','latex')
        xlabel('$x$ (m)','Interpreter','latex')
        ylabel('$y$ (m)','Interpreter','latex')
        zlabel('Amplitude','Interpreter','latex')
        drawnow
        
    end

end

% Get time avarage of the power through each surface
S_time = S_time/i_steps;

%% Simulate in frequency domain

[ebow_freq, hbow_freq] = solve_helmholtz_2d_fd(msh, W, MAT.c, MAT.meps, MAT.mmui, jsbow, -1j * e_exi * E * lz, f, bc);


%% Calculate power emitted through each surface

S_freq = CalcPoyntinvectorXY(msh, ebow_freq, hbow_freq, MAT.ds, MAT.dst);


%% Indices for plots along x-axis

idx_2_plot1 = 1 + ((idx_xmLxh:idx_xpLxh) - 1) * Mx;
idx_2_plot2 = 1 + ((idx_xmLxh:idx_xpLxh) - 1) * Mx + (idx_y0 - 1) * My;


%% Calculate analytical solution for comparision

[E_ana, H_ana] = IncidentPlaneWaveAnalytic(xmesh(idx_2_plot1), -1j * E, f, 1, 4, MAT);


%% Plot analytical solution and solutions in TD and FD

% Plot electric field
figure(2)
subplot(2,1,1);
plot(xmesh(idx_2_plot1), real(E_ana), xmesh(idx_2_plot1), real(ebow_freq(idx_2_plot2+2*np))/lz, xmesh(idx_2_plot1), ebow_old(idx_2_plot2+2*np)/lz);
legend({'Analytical', 'Numerical FD', 'Numerical TD'},'Location','southwest');
title('Electric field in z-direction over x')

% Plot magnetic field
subplot(2,1,2);
plot(xmesh(idx_2_plot1), -real(H_ana), xmesh(idx_2_plot1)-le/2,real(hbow_freq(idx_2_plot2+1*np))/le,xmesh(idx_2_plot1)-le/2,hbow_old(idx_2_plot2+1*np)/le);
legend({'Analytical', 'Numerical FD', 'Numerical TD'},'Location','southwest');
title('Magnetic field in y-direction over x')
drawnow
