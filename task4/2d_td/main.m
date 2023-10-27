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


%% PART 1: General settings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%% Define important parameters for all simulations

% Polarization: 1 for z and 2 for y
polarization = 2;

% Elements per wavelength
elem_per_wavelength = 7;


%% PART 2: Single Slit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define important parameters for the simulation

% Offset in each direction in elements
offset = [0,3,3,3]*4*elem_per_wavelength;
bc.NPML = offset-5;

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
x_offset1 = (-offset(4):-1) * 2 * le;
x_offset2 = L * 1e-6 + (1:offset(2)) * 2 * le;
x_basic = linspace(0, L * 1e-6, points_x);
xmesh = [x_offset1, x_basic, x_offset2];

% Calculate ymesh with respect to the choosen offset
y_offset1 = (-offset(1):-1) * 2 * le;
y_offset2 = h/2 * 1e-6 + (1:offset(3)) * 2 * le;
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
idx_x0      = find(round(xmesh*1e6,4) == 0);

% Index of x = L
idx_xL      = find(round(xmesh*1e6,4) == L);

% Index of y = 0
idx_y0      = find(round(ymesh*1e6,4) == 0);

% Index of y = delta/2
idx_yd_h    = find(round(ymesh*1e6,4) == delta/2);

% Index of y = h/2
idx_h_h     = find(round(ymesh*1e6,4) == h/2);


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


%% Calculate conductivity matrix for conductive PML (open boundary)

[MAT] = conductivePML_2D(bc, msh, MAT, f2);


%% Set up parameters for the simulation in time domain for excitation 1 

% Time step size
dt = CFL(msh, MAT);

% Fit dt to period of excitation 1
dt = 1/f1 / ceil(1/f1 / dt);

% Calculate end time regarding the time needed for the calculation of the
% on avarage emitted power
t_end = sqrt((L * 1e-6)^2 + ((h/2 + delta/2) * 1e-6)^2) / sqrt(material_regions.mu0i/material_regions.epsilon0);
t_end = t_end + 6/f1;


%% Simulate in time domain for excitation 1 

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
    [ebow_new,hbow_new] = solve_FullLeapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.kaps,MAT.c,dt,W);

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
if polarization == 1
    title('Electric field in z-direction (430nm) - Single slit','Interpreter','latex');
elseif polarization == 2
    title('Electric field in y-direction (430nm) - Single slit','Interpreter','latex');
end
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
zlabel('Electric field in V/m','Interpreter','latex')
drawnow


%% Plot avarage power on screen for excitation 1

% Calculate indices of the screen
n_screen = 1 + (idx_xL - 1) * Mx + ((idx_y0:idx_h_h-1) - 1) * My;

% Calculate y-coordinates of corresponding faces as parts of the screen
y_coord_screen_single = ymesh(idx_y0:idx_h_h-1) + le/2;

% Get analytical solution to compare
[I1_ana, bright, dark] = intensityCalcSignleSlit(max(S_ex1(n_screen + np)), delta * 1e-6, L * 1e-6, lambda_1, y_coord_screen_single);

% Vector with numerical calculated intensity
I1_num = S_ex1(n_screen + np);

% Plot avarage power on screen over associated y-coordinates
figure(2)
plot(y_coord_screen_single, I1_num, 'Color', [.5 0 .5], LineWidth=1.5);
hold on;
plot(y_coord_screen_single, I1_ana, 'Color', [0 0 0], LineWidth=1.5, LineStyle='--');
hold on;
scatter(bright, zeros(max(size(bright)),1), "green", 'filled');
hold on;
scatter(dark, zeros(max(size(dark)),1), "black", 'filled');
grid on;
xlabel('$y$ (m)','Interpreter','latex');
ylabel('$Intensity$ (W/(m**2))','Interpreter','latex');
legend({'Numerical solution', 'Analytical solution', 'Loc bright fringes', 'Loc dark fringes'},'Location','northeast');
title('Intensity on screen for wavelength 430nm - Single slit');
drawnow


%% Set up parameters for the simulation in time domain for excitation 2 

% Time step size
dt = CFL(msh, MAT);

% Fit dt to period of excitation 2
dt = 1/f2 / ceil(1/f2 / dt);

% Calculate end time regarding the time needed for the calculation of the
% on avarage emitted power
t_end = sqrt((L * 1e-6)^2 + (h/2 * 1e-6)^2) / sqrt(material_regions.mu0i/material_regions.epsilon0);
t_end = t_end + 6/f2;


%% Simulate in time domain for excitation 2 

% Initialize ebow and hbow
ebow_new = zeros(3*np,1);
hbow_new = zeros(3*np,1);

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
    [ebow_new,hbow_new] = solve_FullLeapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.kaps,MAT.c,dt,W);

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
if polarization == 1
    title('Electric field in z-direction (510nm) - Single slit','Interpreter','latex');
elseif polarization == 2
    title('Electric field in y-direction (510nm) - Single slit','Interpreter','latex');
end
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
zlabel('Electric field in V/m','Interpreter','latex')
drawnow


%% Plot avarage power on screen for excitation 2

% Get analytical solution to compare
[I2_ana, bright, dark] = intensityCalcSignleSlit(max(S_ex2(n_screen + np)), delta * 1e-6, L * 1e-6, lambda_2, y_coord_screen_single);

% Vector with numerically calculated intensity
I2_num = S_ex2(n_screen + np);

% Plot avarage power on screen over associated y-coordinates
figure(4)
plot(y_coord_screen_single, I2_num, 'Color', [0 1 1], LineWidth=1.5);
hold on;
plot(y_coord_screen_single, I2_ana, 'Color', [0 0 0], LineWidth=1.5, LineStyle='--');
scatter(bright, zeros(max(size(bright)),1), "green", 'filled');
hold on;
scatter(dark, zeros(max(size(dark)),1), "black", 'filled');
grid on;
xlabel('$y$ (m)','Interpreter','latex');
ylabel('$Intensity$ (W/(m**2))','Interpreter','latex');
legend({'Numerical solution', 'Analytical solution', 'Loc bright fringes', 'Loc dark fringes'},'Location','northeast');
title('Intensity on screen for wavelength 510nm - Single slit');
drawnow


%% Calculate errors of the solutions for both wavelengths
error_I1 = norm(I1_num - I1_ana') / norm(I1_ana);
error_I2 = norm(I2_num - I2_ana') / norm(I2_ana);

disp(['Relative error of intensity for wavelength 430nm: ', num2str(error_I1)]);
disp(['Relative error of intensity for wavelength 510nm: ', num2str(error_I2)]);


%% Plot the intensity of both wavelengths together
I_combine = I1_num + I2_num;
figure(5)
plot(y_coord_screen_single, I1_num, 'Color', [.5 0 .5], 'LineWidth', 2);
hold on;
plot(y_coord_screen_single, I2_num, 'Color', [0 1 1], 'LineWidth', 2);
hold on;
plot(y_coord_screen_single, I_combine, 'Color', [0 0 0], 'LineStyle','--', LineWidth=1.5);
grid on;
xlabel('$y$ (m)','Interpreter','latex');
ylabel('$Intensity$ (W/(m**2))','Interpreter','latex');
legend({'Wave 430nm', 'Wave 510nm', 'Combined intensity'},'Location','northeast');
title('Combinded intensity of both light waves (430nm and 510nm) - Single slit');
drawnow


%% PART 3: Double slit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clear old variables

clear MAT  material_regions offset bc boxesKappa;


%% Define important parameters for the simulation

% Offset in each direction in wavelengths
offset = [1,1,1,1]*4*elem_per_wavelength+5;

% Edit boundary conditions
bc.bc = ["OPEN", "OPEN", "OPEN", "OPEN"];
bc.NPML = offset-5;


%% Edit basic calculation domain 

% Distance in x-direction in micro meters
L = 10;

% Screen height in micro meters
h = 40;

% Width of slit in micro meters
delta = 1;

% Distance between the slits in micro meters
dist = 3;

% Define maximal edge lenth in the mesh
le = min(lambda_1,lambda_2)/elem_per_wavelength;

% Fit edge length to be an integer divisor of 0.5 micro meter and calculate
% number of edges per micro meter
num_e   = 2 * ceil(0.5e-6 / le);
le      = 1e-6 / num_e;

% Calculate number of points in x- and y-direction
points_x = num_e * L + 1;
points_y = num_e * h + 1;

% Calculate xmesh with respect to the choosen offset
x_offset1 = (-offset(4):-1) * 2 * le;
x_offset2 = L * 1e-6 + (1:offset(2)) * 2 * le;
x_basic = linspace(0, L * 1e-6, points_x);
xmesh = [x_offset1, x_basic, x_offset2];

% Calculate ymesh with respect to the choosen offset
y_offset1 = -h/2 * 1e-6 + (-offset(1):-1) * 2 * le;
y_offset2 = h/2 * 1e-6 + (1:offset(3)) * 2 * le;
y_basic = linspace(-h/2 * 1e-6, h/2 * 1e-6, points_y);
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
idx_x0      = find(round(xmesh*1e6,4) == 0);

% Index of x = L
idx_xL      = find(round(xmesh*1e6,4) == L);

% Index of y = 0
idx_y0      = find(round(ymesh*1e6,4) == 0);

% Index of y = -dist/2 
idx_ymdisth = find(round(ymesh*1e6,4) == -dist/2);

% Index of y = dist/2 
idx_ypdisth = find(round(ymesh*1e6,4) == dist/2);

% Index of y = -dist/2 - delta
idx_ymdd    = find(round(ymesh*1e6,4) == -dist/2 - delta);

% Index of y = dist/2 + delta
idx_ypdd    = find(round(ymesh*1e6,4) == dist/2 + delta);

% Index of y = -h/2
idx_y_mhh   = find(round(ymesh*1e6,4) == -h/2);

% Index of y = h/2
idx_y_phh   = find(round(ymesh*1e6,4) == h/2);


%% Create the vectors describing the excitation

% Create empty jsbow_space vector
jsbow_excitation = NaN(3*np, 1);

% Create excitation vector for the electric field
e_exitation = NaN(3*np, 1);

% Generate excitation for polarization in z-direction
if polarization == 1

    % Determine indices for points in the double slit
    n_exi1 = 1 + (idx_x0 - 1) * Mx + ((idx_ymdd:idx_ymdisth) - 1) * My;
    n_exi2 = 1 + (idx_x0 - 1) * Mx + ((idx_ypdisth:idx_ypdd) - 1) * My;

    % Use corresponding edges in z-direction for excitation
    e_exitation([n_exi1, n_exi2] + 2*np) = lz;

% Generate excitation for polarization in y-direction
elseif polarization == 2

    % Determine indices for points in the double slit
    n_exi1 = 1 + (idx_x0 - 1) * Mx + ((idx_ymdd:idx_ymdisth-1) - 1) * My;
    n_exi2 = 1 + (idx_x0 - 1) * Mx + ((idx_ypdisth:idx_ypdd-1) - 1) * My;

    % Use corresponding edges in y-direction for excitation
    e_exitation([n_exi1, n_exi2] + 1*np) = le; 

end


%% Edit material regions

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


%% Apply boundary conditions and get excitation vectors for the simulation

[bc, W, e_exi, jsbow] = apply_bc(msh, bc, e_exitation, jsbow_excitation);


%% Generate matrices for calculation

[MAT, bc] = generate_MAT(msh, bc, material_regions, ["CurlP"]); %#ok<NBRAK2> 


%% Apply open boundary condition with conducting PML boundary

MAT = conductivePML_2D(bc, msh, MAT, f2);


%% Set up parameters for the simulation in time domain for excitation 1

% Time step size
dt = CFL(msh, MAT);

% Fit dt to period of excitation 1
dt = 1/f1 / ceil(1/f1 / dt);

% Calculate end time regarding the time needed for the calculation of the
% on avarage emitted power
t_end = sqrt((L * 1e-6)^2 + (h/2 * 1e-6 + dist/2 * 1e-6 + delta * 1e-6)^2) / sqrt(material_regions.mu0i/material_regions.epsilon0);
t_end = t_end + 6/f1;


%% Simulate in time domain for excitation 1 

% Initialize ebow and hbow
ebow_new = zeros(3*np,1);
hbow_new = zeros(3*np,1);

% Add inverse permittivity matrix
MAT.mepsi = nullInv(MAT.meps);

% Initialize calculation of avarage power
S_ex1d = zeros(3*np,1);
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
    [ebow_new,hbow_new] = solve_FullLeapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.kaps,MAT.c,dt,W);

    % Sum up power through each surface for one period before the end
    if t >= t_end - 1/f1
        S_ex1d = S_ex1d + CalcPoyntinvectorXY(msh, ebow_new, hbow_new, MAT.ds, MAT.dst);
        i_steps = i_steps + 1;
    end

end

% Get time avarage of the power through each surface
S_ex1d = S_ex1d/i_steps;


%% Plot electric field over domain for excitation 1

figure(6)
[X,Y] = meshgrid(xmesh, ymesh);
if polarization == 1
    e_surf = reshape(ebow_new(2*np+1:3*np,1)/lz, [msh.nx, msh.ny]);
elseif polarization == 2
    e_surf = reshape(ebow_new(np+1:2*np,1)/le, [msh.nx, msh.ny]);
end
e_surf_plot = surf(X,Y,e_surf');
xlim([0, L*1e-6])
ylim([-h/2*1e-6, h/2*1e-6])
set(e_surf_plot,'LineStyle','none')
view(2)
colormap winter;
if polarization == 1
    title('Electric field in z-direction (430nm) - Double slit','Interpreter','latex');
elseif polarization == 2
    title('Electric field in y-direction (430nm) - Double slit','Interpreter','latex');
end
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
zlabel('Electric field in V/m','Interpreter','latex')
drawnow


%% Plot avarage power on screen for excitation 1

% Calculate indices of the screen
n_screen = 1 + (idx_xL - 1) * Mx + ((idx_y_mhh:idx_y_phh-1) - 1) * My;

% Calculate y-coordinates of corresponding faces as parts of the screen
y_coord_screen_double = ymesh(idx_y_mhh:idx_y_phh-1) + le/2;

% Vector with numerical calculated intensity
I1d_num = S_ex1d(n_screen + np);

% Calculate intensity spectrum to correspondinng single slit
I1_single = [flip(I1_num) ; I1_num(2:end,1)];
y1_single = [-flip(y_coord_screen_single), y_coord_screen_single(1, 2:end)];

% Calculate intensity spectrum without interference
I11 = interp1(y1_single - (dist/2 + delta/2) * 1e-6, I1_single, y_coord_screen_double, 'pchip', 0);
I12 = interp1(y1_single + (dist/2 + delta/2) * 1e-6, I1_single, y_coord_screen_double, 'pchip', 0);
I1_double_single = I11 + I12 + 2 * I11.^0.5 .* I12.^0.5;

% Plot avarage power on screen over associated y-coordinates
figure(7)
plot(y_coord_screen_double, I1d_num, 'Color', [.5 0 .5], LineWidth=1.5);
hold on;
plot(y_coord_screen_double, I1_double_single, 'Color', [0 0 0], 'LineStyle','--', LineWidth=1.5);
xlabel('$y$ (m)','Interpreter','latex');
ylabel('$Intensity$ (W/(m**2))','Interpreter','latex');
legend({'Numerical solution', 'Envelope two single slits'},'Location','northeast');
title('Intensity on screen for wavelength 430nm - Double slit');
drawnow


%% Set up parameters for the simulation in time domain for excitation 2

% Time step size
dt = CFL(msh, MAT);

% Fit dt to period of excitation 2
dt = 1/f2 / ceil(1/f2 / dt);

% Calculate end time regarding the time needed for the calculation of the
% on avarage emitted power
t_end = sqrt((L * 1e-6)^2 + (h/2 * 1e-6 + dist/2 * 1e-6 + delta * 1e-6)^2) / sqrt(material_regions.mu0i/material_regions.epsilon0);
t_end = t_end + 6/f2;


%% Simulate in time domain for excitation 2

% Initialize ebow and hbow
ebow_new = zeros(3*np,1);
hbow_new = zeros(3*np,1);

% Initialize calculation of avarage power
S_ex2d = zeros(3*np,1);
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
    [ebow_new,hbow_new] = solve_FullLeapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.kaps,MAT.c,dt,W);

    % Sum up power through each surface for one period before the end
    if t >= t_end - 1/f2
        S_ex2d = S_ex2d + CalcPoyntinvectorXY(msh, ebow_new, hbow_new, MAT.ds, MAT.dst);
        i_steps = i_steps + 1;
    end

end

% Get time avarage of the power through each surface
S_ex2d = S_ex2d/i_steps;


%% Plot electric field over domain for excitation 2

figure(8)
[X,Y] = meshgrid(xmesh, ymesh);
if polarization == 1
    e_surf = reshape(ebow_new(2*np+1:3*np,1)/lz, [msh.nx, msh.ny]);
elseif polarization == 2
    e_surf = reshape(ebow_new(np+1:2*np,1)/le, [msh.nx, msh.ny]);
end
e_surf_plot = surf(X,Y,e_surf');
xlim([0, L*1e-6])
ylim([-h/2*1e-6, h/2*1e-6])
set(e_surf_plot,'LineStyle','none')
view(2)
colormap winter;
if polarization == 1
    title('Electric field in z-direction (510nm) - Double slit','Interpreter','latex');
elseif polarization == 2
    title('Electric field in y-direction (510nm) - Double slit','Interpreter','latex');
end
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
zlabel('Electric field in V/m','Interpreter','latex')
drawnow


%% Plot avarage power on screen for excitation 2

% Vector with numerical calculated intensity
I2d_num = S_ex2d(n_screen + np);

% Calculate intensity spectrum to corresponding single slit
I2_single = [flip(I2_num) ; I2_num(2:end,1)];
y2_single = [-flip(y_coord_screen_single), y_coord_screen_single(1, 2:end)];

% Calculate intensity spectrum without interference
I21 = interp1(y2_single - (dist/2 + delta/2) * 1e-6, I2_single, y_coord_screen_double, 'pchip', 0);
I22 = interp1(y2_single + (dist/2 + delta/2) * 1e-6, I2_single, y_coord_screen_double, 'pchip', 0);
I2_double_single = I21 + I22 + 2 * I21.^0.5 .* I22.^0.5;

% Plot avarage power on screen over associated y-coordinates
figure(9)
plot(y_coord_screen_double, I2d_num, 'Color', [0 1 1], LineWidth=1.5);
hold on;
plot(y_coord_screen_double, I2_double_single, 'Color', [0 0 0], 'LineStyle','--', LineWidth=1.5);
hold on;
xlabel('$y$ (m)','Interpreter','latex');
ylabel('$Intensity$ (W/(m**2))','Interpreter','latex');
legend({'Numerical solution', 'Envelope two single slits'},'Location','northeast');
title('Intensity on screen for wavelength 510nm - Double slit');
drawnow


%% Plot the intensity of both wavelengths together
Id_combine = I1d_num + I2d_num;
envelope_combine = I1_double_single + I2_double_single;
figure(10)
plot(y_coord_screen_double, I1d_num, 'Color', [.5 0 .5], 'LineWidth', 2);
hold on;
plot(y_coord_screen_double, I2d_num, 'Color', [0 1 1], 'LineWidth', 2);
hold on;
plot(y_coord_screen_double, Id_combine, 'Color', [0 0 0], 'LineStyle','--', LineWidth=1.5);
grid on;
plot(y_coord_screen_double, envelope_combine, 'Color', [1 0 0], 'LineStyle','--', LineWidth=1);
xlabel('$y$ (m)','Interpreter','latex');
ylabel('$Intensity$ (W/(m**2))','Interpreter','latex');
legend({'Wave 430nm', 'Wave 510nm', 'Combined intensity', 'Envelope'},'Location','northeast');
title('Combinded intensity of both light waves (430nm and 510nm) - Double slit');
drawnow
