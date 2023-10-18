%% Verification script

% Polarization: Y           Propagation: X

% Verification script for the simulation of a plane wave with polarization
% of the electric field in y-direction and propagation in x-direction. The
% domain is devided in two region. For x<0, the relative permittivity is 1,
% for x>0 it is 4. The simulation is executet in time domain 
% and the results are compared with an analytic solution.

% The domain is a 4 micro meter (x) times 2 micro meter (y) rectangle (2D) 
% and the plane wave for the exitation is placed at x = -3 micro meter. 
% From there, the wave travels to a screen at x = 1 micro meter 
% (no refelction at the screen), where the intensity of the incoming wave 
% is observed. The wave has a wavelength of 430nm and an amplitude of 500V. 


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

% Calculate value of excitation for t=0 at -3 micro meters
k1 = 2*pi*f * sqrt(material_regions.epsilon0/material_regions.mu0i);
Z1 = sqrt(1/material_regions.mu0i/material_regions.epsilon0);
E0 = E * exp(-1j * k1 * (-3e-6));
phi1 = angle(E0); % -> Start angle for excitation in time domain

% Function for excitation in time domain
func_exi    = @(t)(E * cos(2*pi*f*t + phi1));


%% Define important parameters for the simulation

% Elements per wavelength
elem_per_wavelength = 30;

% Offset in each direction in elements
offset = [0,25,0,25];

% Edit boundary conditions
bc.bc = ["PEC", "OPEN", "PEC", "OPEN"];


%% Edit basic calculation domain 

% Distance in x-direction in micro meters
Lx = 4;

% Distance in y-direction in micro meters
Ly = 2;

% Define maximal edge lenth in the mesh
le = lambda/elem_per_wavelength;

% Fit edge length to be an integer divisor of 1 micro meter and calculate
% number of edges per micro meter
num_e   = ceil(1e-6 / le);
le      = 1e-6 / num_e;

% Calculate xmesh with respect to the choosen offset
x_offset1 = -Lx*3/4 * 1e-6  + (-offset(4):-1) * le;
x_offset2 = Lx*1/4 * 1e-6 + (1:offset(2)) * le/2;
x_basic1 = linspace(-Lx*3/4 * 1e-6, 0 - le, 3*num_e);
x_basic2 = linspace(0, Lx*1/4 * 1e-6, 2*num_e+1);
xmesh = [x_offset1, x_basic1, x_basic2, x_offset2];

% Calculate ymesh with respect to the choosen offset
y_offset1 = -Ly/2 * 1e-6  + (-offset(1):-1) * le;
y_offset2 = Ly/2 * 1e-6 + (1:offset(3)) * le;
y_basic = linspace(-Ly/2 * 1e-6, Ly/2 * 1e-6, 2*num_e+1);
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

% Index of x = -Lx*3/4
idx_xmin    = find(xmesh == -Lx*3/4 * 1e-6);

% Index of x = Lx*1/4
idx_xmax    = find(xmesh == Lx*1/4 * 1e-6);

% Index of y = 0
idx_y0      = find(ymesh == 0);

% Index of y = -Ly/2
idx_ymin    = find(ymesh == -Ly/2 * 1e-6);

% Index of y = Ly/2
idx_ymax    = find(ymesh == Ly/2 * 1e-6);


%% Create the vectors describing the excitation

% Create empty jsbow_space vector
jsbow_excitation = NaN(3*np, 1);

% Create excitation vector for the electric field
e_exitation = NaN(3*np, 1);

% Determine indices for points at x = -Lx/2
n_exi = 1 + (idx_xmin - 1) * Mx + ((idx_ymin:idx_ymax-1) - 1) * My;

% Use corresponding edges in y-direction for excitation
e_exitation(n_exi + 1*np) = 1;


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


%% Set up parameters for the simulation in time domain for excitation 1 

% Time step size
dt = 1 * CFL(msh, MAT);

% Fit dt to period of excitation 1
dt = 1/f / ceil(1/f / dt);

% End time
t_end = 14/f;


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
for t = 0:dt:t_end

    % Old values
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    % Calculate value for excitation
    e_exi_old = e_exi * func_exi(t) * le;
    e_exi_new = e_exi * func_exi(t+1*dt) * le;

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

    if steps_movie == 20

        steps_movie = 0;

        [X,Y] = meshgrid(xmesh, ymesh);
        e_surf = reshape(ebow_new(1*np+1:2*np,1)/le, [msh.nx, msh.ny]);
        
        figure(1)
        e_surf_plot = surf(X,Y,e_surf');
    
        xlim([-Lx*3/4*1e-6, Lx*1/4*1e-6])
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


%% Indices for plots along x-axis

idx_2_plot1 = 1 + ((idx_xmin:idx_xmax) - 1) * Mx;
idx_2_plot2 = 1 + ((idx_xmin:idx_xmax) - 1) * Mx + (idx_y0 - 1) * My;


%% Indices for plot on screen

idx_screen_1 = idx_ymin:idx_ymax-1;
idx_screen_2 = 1 + (idx_xmax - 1) * Mx + ((idx_ymin:idx_ymax-1) - 1) * My;


%% Calculate analytical solution for comparision

% Calculate field values
[E_ana, ~] = Analytic_Pol_y_Prop_x(xmesh(idx_2_plot1), E, f, 1, 4, MAT);
[~, H_ana] = Analytic_Pol_y_Prop_x(xmesh(idx_2_plot1) + [diff(xmesh(idx_2_plot1))/2, le/4], E, f, 1, 4, MAT);

% Calculate Poyntinvector
S_ana = real(0.5 * E_ana(end) * conj(H_ana(end)));


%% Plot analytical solution and solutions in TD and FD

% Plot electric field
figure(2)
subplot(2,1,1);
plot(xmesh(idx_2_plot1), real(E_ana), xmesh(idx_2_plot1), ebow_old(idx_2_plot2+1*np)/le);
xlim([-Lx*3/4*1e-6, Lx*1/4*1e-6]);
ylim([-550, 550]);
xlabel('$x$ (m)','Interpreter','latex');
ylabel('Electric field strength','Interpreter','latex');
legend({'Analytical', 'Numerical TD'},'Location','southwest');
title('Electric field in z-direction over x');

% Plot magnetic field
subplot(2,1,2);
plot(xmesh(idx_2_plot1) + [diff(xmesh(idx_2_plot1))/2, le/4], real(H_ana), xmesh(idx_2_plot1) + [diff(xmesh(idx_2_plot1))/2, le/4], hbow_old(idx_2_plot2+2*np)/lz);
xlim([-Lx*3/4*1e-6, Lx*1/4*1e-6]);
ylim([-800/Z1, 800/Z1]);
xlabel('$x$ (m)','Interpreter','latex');
ylabel('Magnetic field strength','Interpreter','latex');
legend({'Analytical', 'Numerical TD'},'Location','southwest');
title('Magnetic field in y-direction over x');


%% Plot solution for power emitted on screen at x = 1 micro meter
figure(3)
plot(ymesh(idx_screen_1)+le/2, ones(max(size(idx_screen_1)),1) * S_ana, ymesh(idx_screen_1)+le/2, S_time(idx_screen_2+np))
xlim([-1*1e-6, 1*1e-6])
ylim([200, 350])
xlabel('$y$ (m)','Interpreter','latex');
ylabel('Avarage power in w/m**2 ','Interpreter','latex');
legend({'Analytical', 'Numerical TD'},'Location','southwest');
title('Avarage power emitted on screen at x = 1 micro meter over y');
drawnow


%% Calculate and display errors

% Errors electric field
error_E_TD = norm(ebow_old(idx_2_plot2+1*np)/le - real(E_ana)') / norm(real(E_ana));

% Errors magnetic field
error_H_TD = norm(hbow_old(idx_2_plot2+2*np)/lz- real(H_ana)') / norm(real(H_ana));

% Error power on screen
error_S_TD = norm(S_time(idx_screen_2+np)- (ones(max(size(idx_screen_1)),1) * S_ana)') / norm(ones(max(size(idx_screen_1)),1) * S_ana);

% Display errors:
disp(['Relative error of the electric field: ', num2str(error_E_TD)]);
disp(['Relative error of the magnetic field: ', num2str(error_H_TD)]);
disp(['Relative error of the emitted power: ', num2str(error_S_TD)]);

