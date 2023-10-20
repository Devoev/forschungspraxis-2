%% Add paths

% Clear variables
clc
clear 
close all

% Get parent directory
filePath = matlab.desktop.editor.getActiveFilename;
[parent] = fileparts(filePath);
parent = fileparts(parent) ;
parent = fileparts(parent);

% Paths to add
path_msh_func = append(parent, '\fit\2d\mesh');
path_mat_func = append(parent, '\fit\2d\matrices');
path_solver_func = append(parent, '\fit\2d\solver');
path_util_func = append(parent, '\fit\2d\util');

% Add paths
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func)


%% Define excitation

% Add basic constants to material_regions object
material_regions.epsilon0 = 8.854187e-12;
material_regions.mu0i = 1/(pi*4e-7);

% Excitation 1
lambda_1    = 430e-9;
f1          = sqrt(material_regions.mu0i/material_regions.epsilon0)/lambda_1;           
E1          = 250;
func_exi_1  = @(t)(E1 * sin(2*pi*f1*t));

% Excitation 2
lambda_2    = 510e-9;
f2          = sqrt(material_regions.mu0i/material_regions.epsilon0)/lambda_2;           
E2          = 500;
func_exi_2  = @(t)(E2 * sin(2*pi*f2*t));


%% Define important parameters for the simulation

% Polarization: 1 for z and 2 for y
polarization = 1;

% Elements per wavelength
elem_per_wavelength = 8;

% Offset in each direction in elements
offset = [0,10,10,10];

% Edit boundary conditions
if polarization == 1
    bc.bc = ["PMC", "OPEN", "OPEN", "OPEN"];
elseif polarization == 2
    bc.bc = ["PEC", "OPEN", "OPEN", "OPEN"];
end


%% Edit basic calculation domain 

% Distance in x-direction
L = 10;

% Screen height
h = 40;

% Width of slit
delta = 1;

% Define maximal edge lenth in the mesh
le = min(lambda_1,lambda_2)/elem_per_wavelength;

% Fit edge length to be an integer divisor of 0.5 micro meter and calculate
% number of edges per micro meter
num_e   = 2 * ceil(0.5e-6 / le);
le      = 1e-6 / num_e;

% Calculate number of points in x- and y-direction
points_x = num_e * L + 1;
points_y = num_e * h/2 + 1;

% Calculate xmesh with respect to the choosen offset
x_offset1 = (-offset(4):-1) * le;
x_offset2 = L * 1e-6 + (1:offset(2)) * le;
x_basic = linspace(0, L * 1e-6, points_x);
xmesh = [x_offset1, x_basic, x_offset2];

% Calculate ymesh with respect to the choosen offset
y_offset1 = (-offset(1):-1) * le;
y_offset2 = h/2 * 1e-6 + (1:offset(3)) * le;
y_basic = linspace(0, h/2 * 1e-6, points_y);
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

% Index of x = L
idx_xL      = find(xmesh == L * 1e-6);

% Index of y = 0
idx_y0      = find(ymesh == 0);

% Index of y = delta/2
idx_yd_h    = find(ymesh == delta/2 * 1e-6);

% Index of y = h/2
idx_h_h     = find(ymesh == h/2 * 1e-6);


%% Create the vectors describing the excitation

% Create empty jsbow_space vector
jsbow_excitation = NaN(3*np, 1);

% Create excitation vector for the electric field
e_exitation = NaN(3*np, 1);

% Generate excitation for polarization in z-direction
if polarization == 1

    % Determine indices for points in the single slit
    n_exi = 1 + (idx_x0 - 1) * Mx + ((idx_y0:idx_yd_h) - 1) * My;

    % Use corresponding edges in z-direction for excitation
    e_exitation(n_exi + 2*np) = lz;

% Generate excitation for polarization in y-direction
elseif polarization == 2

    % Determine indices for points in the single slit
    n_exi = 1 + (idx_x0 - 1) * Mx + ((idx_y0:idx_yd_h-1) - 1) * My;

    % Use corresponding edges in y-direction for excitation
    e_exitation(n_exi + 1*np) = le; 

end


%% Edit material regions and add them to the object material_regions

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


%% Apply boundary conditions and get excitation vectors for the simulation

[bc, W, e_exi, jsbow] = apply_bc(msh, bc, e_exitation, jsbow_excitation);


%% Generate matrices for calculation

[MAT, bc] = generate_MAT(msh, bc, material_regions, ["CurlP"]); %#ok<NBRAK2> 


%% Set up parameters for the simulation in time domain for excitation 1 

% Time step size
dt = CFL(msh, MAT);

% Fit dt to period of excitation 1
dt = 1/f1 / ceil(1/f1 / dt);

% Calculate end time regarding the time needed for the calculation of the
% on avarage emitted power
t_end = sqrt((L * 1e-6)^2 + (h/2 * 1e-6)^2) / sqrt(material_regions.mu0i/material_regions.epsilon0);
t_end = t_end + 3/f1;


%% Simulate in time domain for excitation 1 

% Initialize open boundary condition if needed
[mur_edges,mur_n_edges, mur_deltas] = initMur_2D(msh, bc);

% Initialize ebow and hbow
ebow_new = zeros(3*np,1);
hbow_new = zeros(3*np,1);

% Add inverse permittivity matrix
MAT.mepsi = nullInv(MAT.meps);

% Initialize calculation of avarage power
S_ex1 = zeros(3*np,1);
i_steps = 0;

% Calculate time steps
for t = 0:dt:t_end

    % Old values
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    % Calculate value for excitation
    e_exi_old = e_exi * func_exi_1(t);
    e_exi_new = e_exi * func_exi_1(t+1);

    % Execute timestep with leapfrog
    [ebow_new,hbow_new] = solve_leapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.c,dt,W);

    % Apply open boundary with mur condition
    ebow_new = applyMur_2D(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt, bc);

    % Sum up power through each surface for one period before the end
    if t >= t_end - 1/f1
        S_ex1 = S_ex1 + CalcPoyntinvectorXY(msh, ebow_new, hbow_new, MAT.ds, MAT.dst);
        i_steps = i_steps + 1;
    end

end

% Get time avarage of the power through each surface
S_ex1 = S_ex1/i_steps;


%% Plot electric field over domain for excitation 1

figure(1)
[X,Y] = meshgrid(xmesh, ymesh);
if polarization == 1
    e_surf = reshape(ebow_new(2*np+1:3*np,1)/lz, [msh.nx, msh.ny]);
elseif polarization == 2
    e_surf = reshape(ebow_new(np+1:2*np,1)/le, [msh.nx, msh.ny]);
end
e_surf_plot = surf(X,Y,e_surf');
xlim([0, L*1e-6])
ylim([0, h/2*1e-6])
set(e_surf_plot,'LineStyle','none')
view(2)
colormap winter;
title('Value of electric field in z-direction (430nm)','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
zlabel('Electric field in V/m','Interpreter','latex')
drawnow


%% Plot avarage power on screen for excitation 1

% Calculate indices of the screen
n_screen = 1 + (idx_xL - 1) * Mx + ((idx_y0:idx_h_h-1) - 1) * My;

% Calculate y-coordinates of corresponding faces as parts of the screen
y_coord_screen = ymesh(idx_y0:idx_h_h-1) + le/2;

% Get analytical solution to compare
[I1_ana, bright, dark] = intensityCalcSignleSlit(max(S_ex1(n_screen + np)), delta * 1e-6, L * 1e-6, lambda_1, y_coord_screen);

% Vector with numerical calculated intensity
I1_num = S_ex1(n_screen + np);

% Plot avarage power on screen over associated y-coordinates
figure(2)
plot(y_coord_screen, I1_num, 'Color', [.5 0 .5], LineWidth=1.5);
hold on;
plot(y_coord_screen, I1_ana, 'Color', [0 0 0], LineWidth=1.5, LineStyle='--');
hold on;
scatter(bright, zeros(max(size(bright)),1), "green", 'filled');
hold on;
scatter(dark, zeros(max(size(dark)),1), "black", 'filled');
grid on;
xlabel('$y$ (m)','Interpreter','latex');
ylabel('$Intensity$ (W/(m**2))','Interpreter','latex');
legend({'Numerical TD', 'Analytical', 'Loc bright fringes', 'Loc dark fringes'},'Location','northeast');
title('Intensity on screen for electric field with 430nm wavelength');
drawnow


%% Set up parameters for the simulation in time domain for excitation 2 

% Time step size
dt = CFL(msh, MAT);

% Fit dt to period of excitation 2
dt = 1/f2 / ceil(1/f2 / dt);

% Calculate end time regarding the time needed for the calculation of the
% on avarage emitted power
t_end = sqrt((L * 1e-6)^2 + (h/2 * 1e-6)^2) / sqrt(material_regions.mu0i/material_regions.epsilon0);
t_end = t_end + 3/f2;


%% Simulate in time domain for excitation 2 

% Initialize open boundary condition if needed
[mur_edges,mur_n_edges, mur_deltas] = initMur_2D(msh, bc);

% Initialize ebow and hbow
ebow_new = zeros(3*np,1);
hbow_new = zeros(3*np,1);

% Add inverse permittivity matrix
MAT.mepsi = nullInv(MAT.meps);

% Initialize calculation of avarage power
S_ex2 = zeros(3*np,1);
i_steps = 0;

% Calculate time steps
for t = 0:dt:t_end

    % Old values
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    % Calculate value for excitation
    e_exi_old = e_exi * func_exi_2(t);
    e_exi_new = e_exi * func_exi_2(t+1);

    % Execute timestep with leapfrog
    [ebow_new,hbow_new] = solve_leapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.c,dt,W);

    % Apply open boundary with mur condition
    ebow_new = applyMur_2D(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt, bc);

    % Sum up power through each surface for one period before the end
    if t >= t_end - 1/f2
        S_ex2 = S_ex2 + CalcPoyntinvectorXY(msh, ebow_new, hbow_new, MAT.ds, MAT.dst);
        i_steps = i_steps + 1;
    end

end

% Get time avarage of the power through each surface
S_ex2 = S_ex2/i_steps;


%% Plot electric field over domain for excitation 2

figure(3)
[X,Y] = meshgrid(xmesh, ymesh);
if polarization == 1
    e_surf = reshape(ebow_new(2*np+1:3*np,1)/lz, [msh.nx, msh.ny]);
elseif polarization == 2
    e_surf = reshape(ebow_new(np+1:2*np,1)/le, [msh.nx, msh.ny]);
end
e_surf_plot = surf(X,Y,e_surf');
xlim([0, L*1e-6])
ylim([0, h/2*1e-6])
set(e_surf_plot,'LineStyle','none')
view(2)
colormap winter;
title('Value of electric field in z-direction (510nm)','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
zlabel('Electric field in V/m','Interpreter','latex')
drawnow


%% Plot avarage power on screen for excitation 2

% Calculate indices of the screen
n_screen = 1 + (idx_xL - 1) * Mx + ((idx_y0:idx_h_h-1) - 1) * My;

% Calculate y-coordinates of corresponding faces as parts of the screen
y_coord_screen = ymesh(idx_y0:idx_h_h-1) + le/2;

% Get analytical solution to compare
[I2_ana, bright, dark] = intensityCalcSignleSlit(max(S_ex2(n_screen + np)), delta * 1e-6, L * 1e-6, lambda_2, y_coord_screen);

% Vector with numerically calculated intensity
I2_num = S_ex2(n_screen + np);

% Plot avarage power on screen over associated y-coordinates
figure(4)
plot(y_coord_screen, I2_num, 'Color', [0 1 1], LineWidth=1.5);
hold on;
plot(y_coord_screen, I2_ana, 'Color', [0 0 0], LineWidth=1.5, LineStyle='--');
scatter(bright, zeros(max(size(bright)),1), "green", 'filled');
hold on;
scatter(dark, zeros(max(size(dark)),1), "black", 'filled');
grid on;
xlabel('$y$ (m)','Interpreter','latex');
ylabel('$Intensity$ (W/(m**2))','Interpreter','latex');
legend({'Numerical', 'Analytical', 'Loc bright fringes', 'Loc dark fringes'},'Location','northeast');
title('Intensity on screen for electric field with 510nm wavelength');
drawnow


%% Calculate errors of the solutions for both wavelengths
error_I1 = norm(I1_num - I1_ana') / norm(I1_ana);
error_I2 = norm(I2_num - I2_ana') / norm(I2_ana);

disp(['Relative error of intensity for wavelength 430nm: ', num2str(error_I1)]);
disp(['Relative error of intensity for wavelength 510nm: ', num2str(error_I2)]);


%% Plot the intensity of both wavelengths together
I_combine = I1_num + I2_num;
figure(5)
plot(y_coord_screen, I1_num, 'Color', [.5 0 .5], 'LineWidth', 2);
hold on;
plot(y_coord_screen, I2_num, 'Color', [0 1 1], 'LineWidth', 2);
hold on;
plot(y_coord_screen, I_combine, 'Color', [0 0 0], 'LineStyle','--', LineWidth=1.5);
grid on;
xlabel('$y$ (m)','Interpreter','latex');
ylabel('$Intensity$ (W/(m**2))','Interpreter','latex');
legend({'Wave 430nm', 'Wave 510nm', 'Combined intensity'},'Location','northeast');
title('Combinded intensity of both light waves (430nm and 510nm)');
drawnow
