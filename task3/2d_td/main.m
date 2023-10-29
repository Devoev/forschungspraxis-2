%% Task3 - Thin film - time domain solver

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
polarization = 2;

% Elements per wavelength
elem_per_wavelength = 16;

% Offset in each direction
offset = [0,8,0,12]*elem_per_wavelength;

% Edit boundary conditions
if polarization == 1
    bc.bc = ["PMC", "OPEN", "PMC", "OPEN"];
elseif polarization == 2
    bc.bc = ["PEC", "OPEN", "PEC", "OPEN"];
end
bc.NPML = [0,8,0,8]*elem_per_wavelength - 3;


%% Edit basic calculation domain 

% Distance in x-direction in micro meter
L = 10;

% Screen height in micro meter
h = 4;

% Thcikness of thin film in micro meter
a = 0.1;

% Define maximal edge lenth in the mesh
le = min(lambda_1,lambda_2)/elem_per_wavelength;

% Fit edge length to be an integer divisor of 0.1 micro meter and calculate
% number of edges per micro meter
num_e   = 10 * ceil(0.1e-6 / le);
le      = 1e-6 / num_e;

% Calculate number of points in y-direction
points_y = num_e * h + 1;

% Calculate xmesh with respect to the choosen offset
x_offset1 = (-offset(4):-1) * le;
x_offset2 = L * 1e-6 + (1:offset(2)) * le;
x_basic = [linspace(0, (L/2)*1e-6, num_e*(L/2)+1), ...
    linspace((L/2)*1e-6 + le/2, (L/2+a)*1e-6 - le/2, ceil(2*num_e*a)), ...
    linspace((L/2+a)*1e-6, L*1e-6, num_e*(L/2)+1)];
xmesh = [x_offset1, x_basic, x_offset2];

% Calculate ymesh with respect to the choosen offset
y_offset1 = -h/2 * 1e-6 + (-offset(1):-1) * le;
y_offset2 = h/2 * 1e-6 + (1:offset(3)) * le;
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

% Index in x-direction for excitation (inside the offset!!!!)
idx_exi     = bc.NPML(4) + 3;

% Index of x = L/2
idx_xL_h    = find(round(xmesh*1e6,4) == L/2);

% Index of x = L/2 + a
idx_xL_h_a  = find(round(xmesh*1e6,4) == L/2 + a);

% Index of x = L
idx_xL      = find(round(xmesh*1e6,4) == L);

% Index of y = -h/2
idx_ymh_h   = find(round(ymesh*1e6,4) == -h/2);

% Index of y = 0
idx_y0      = find(round(ymesh*1e6,4) == 0);

% Index of y = h/2
idx_yph_h   = find(round(ymesh*1e6,4) == h/2);


%% Create the vectors describing the excitation

% Create empty jsbow_space vector
jsbow_excitation = NaN(3*np, 1);

% Create excitation vector for the electric field
e_exitation = NaN(3*np, 1);

% Generate excitation for polarization in z-direction
if polarization == 1

    % Determine indices for points in the single slit
    n_exi = 1 + (idx_exi - 1) * Mx + ((idx_ymh_h:idx_yph_h) - 1) * My;

    % Use corresponding edges in z-direction for excitation
    e_exitation(n_exi + 2*np) = lz;

% Generate excitation for polarization in y-direction
elseif polarization == 2

    % Determine indices for points in the single slit
    n_exi = 1 + (idx_exi - 1) * Mx + ((idx_ymh_h:idx_yph_h-1) - 1) * My;

    % Use corresponding edges in y-direction for excitation
    e_exitation(n_exi + 1*np) = le; 

end


%% Edit material regions and add them to the object material_regions

% Regions for relative permittivity
boxesEpsilonR(1).box = [1, idx_xL_h, 1, ny];
boxesEpsilonR(1).value = 1;
boxesEpsilonR(2).box = [idx_xL_h, idx_xL_h_a, 1, ny];
boxesEpsilonR(2).value = 4;
boxesEpsilonR(3).box = [idx_xL_h_a, nx, 1, ny];
boxesEpsilonR(3).value = 1;
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

[MAT] = conductivePML_2D(bc, msh, MAT, f1);


%% Set up parameters for the simulation in time domain for excitation 1 

% Time step size
dt = CFL(msh, MAT);

% Fit dt to period of excitation 1
dt = 1/f1 / ceil(1/f1 / dt);

% Calculate end time regarding the time needed for the calculation of the
% on avarage emitted power
% time of wave arrival + four percent (for possible num. errors)
c0 = sqrt(MAT.mu0i/MAT.epsilon0);
t_end = (1/3 * offset(4) *le + L/2 * 1e-6) * 2 / c0 + 3/f1;


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
    e_exi_new = e_exi * func_exi_1(t+dt);

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


%% Plot avarage power on screen at x=0 and x=L

% Calculate indices of screen at x = 0 (y-surfaces)
n_screen_0 = 1 + (idx_x0 - 1) * Mx + ((idx_ymh_h:idx_yph_h-1)-1) * My + msh.np;

% Calculate indices of the second screen at x = L (y-surfaces)
n_screen_L = 1 + (idx_xL - 1) * Mx + ((idx_ymh_h:idx_yph_h-1)-1) * My +msh.np;

% Calculate y-coordinates of corresponding faces as parts of the screen
y_coord_screen = msh.ymesh(idx_ymh_h:idx_yph_h-1) + le/2;

% Calculate analytical solution for excitation 1
[S1_f1, S3_f1] = AnaSolPoyntin(E1, f1,material_regions.epsilon0, 4 * material_regions.epsilon0, material_regions.mu0i, material_regions.mu0i, a*1e-6);
S1_f1 = ones(max(size(y_coord_screen)), 1) * S1_f1;
S3_f1 = ones(max(size(y_coord_screen)), 1) * S3_f1;

% Plot avarage power on screen at x = 0
figure(2)
plot(y_coord_screen, S_ex1(n_screen_0),'DisplayName', 'Numerical', 'color', '#1e8080')
hold on
plot(y_coord_screen, S1_f1, 'DisplayName', 'Analytic', 'color', 'red')
hold on
title('Intensity at $x=0$m for excitation 1','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
legend()
drawnow

% Plot avarage power on screen at x = L
figure(3)
plot(y_coord_screen, S_ex1(n_screen_L), 'DisplayName', 'Numerical', 'color', '#1e8080')
hold on
plot(y_coord_screen, S3_f1, 'DisplayName', 'Analytic', 'color', 'red')
hold on
title('Intensity at $x=L=10^{-6}$m for excitation 1','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
legend()
drawnow

% Plot intensity distribution on domain
figure(4)
[X,Y] = meshgrid(msh.xmesh, msh.ymesh);
s_surf = reshape(S_ex1((1+np:2*np),1), [msh.nx, msh.ny]);
s_surf_plot = surf(X,Y,real(s_surf'));
title('Numeric intensity across domain for excitation 1','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
xlim([min(xmesh), max(xmesh)])
ylim([0, h/2*1e-6 - le])
set(s_surf_plot,'LineStyle','none')
colormap winter;
drawnow


%% Calculate conductivity matrix for conductive PML (open boundary)

[MAT] = conductivePML_2D(bc, msh, MAT, f2);


%% Set up parameters for the simulation in time domain for excitation 2 

% Time step size
dt = CFL(msh, MAT);

% Fit dt to period of excitation 2
dt = 1/f2 / ceil(1/f2 / dt);

% Calculate end time regarding the time needed for the calculation of the
% on avarage emitted power
% time of wave arrival + four percent (for possible num. errors)
c0 = sqrt(MAT.mu0i/MAT.epsilon0);
t_end = (1/3 * offset(4) *le + L/2 * 1e-6) * 2 / c0 + 3/f1;


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
    e_exi_new = e_exi * func_exi_2(t+dt);

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


%% Plot avarage power on screen at x=0 and x=L

% Calculate analytical solution for excitation 2
[S1_f2, S3_f2] = AnaSolPoyntin(E2,f2,material_regions.epsilon0, 4 * material_regions.epsilon0, material_regions.mu0i, material_regions.mu0i, a*1e-6);
S1_f2 = ones(max(size(y_coord_screen)), 1) * S1_f2;
S3_f2 = ones(max(size(y_coord_screen)), 1) * S3_f2;

% Plot avarage power on screen at x = 0
figure(5)
plot(y_coord_screen, S_ex2(n_screen_0),'DisplayName', 'Numerical', 'color', '#1e8080')
hold on
plot(y_coord_screen, S1_f2, 'DisplayName', 'Analytic', 'color', 'red')
hold on
title('Intensity at $x=0$m for excitation 2','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
legend()
drawnow

% Plot avarage power on screen at x = L
figure(6)
plot(y_coord_screen, S_ex2(n_screen_L), 'DisplayName', 'Numerical', 'color', '#1e8080')
hold on
plot(y_coord_screen, S3_f2, 'DisplayName', 'Analytic', 'color', 'red')
hold on
title('Intensity at $x=L=10^{-6}$m for excitation 2','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
legend()
drawnow

% Plot intensity distribution on domain
figure(7)
[X,Y] = meshgrid(msh.xmesh, msh.ymesh);
s_surf = reshape(S_ex2((1+np:2*np),1), [msh.nx, msh.ny]);
s_surf_plot = surf(X,Y,real(s_surf'));
title('Numeric intensity across domain for excitation 2','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
xlim([min(xmesh), max(xmesh)])
ylim([0, h/2*1e-6 - le])
set(s_surf_plot,'LineStyle','none')
colormap winter;
drawnow


%% Plot combined intensities on screen at x=0 and x=L

% Plot avarage power on screen at x = 0
figure(8)
plot(y_coord_screen, S_ex1(n_screen_0) + S_ex2(n_screen_0),'DisplayName', 'Numerical', 'color', '#1e8080')
hold on
plot(y_coord_screen, S1_f1 + S1_f2, 'DisplayName', 'Analytic', 'color', 'red')
hold on
title('Combined intensity at $x=0$m for excitation 1 and 2','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
legend()
drawnow

% Plot avarage power on screen at x = L
figure(9)
plot(y_coord_screen, S_ex1(n_screen_L) + S_ex2(n_screen_L), 'DisplayName', 'Numerical', 'color', '#1e8080')
hold on
plot(y_coord_screen, S3_f1 + S3_f2, 'DisplayName', 'Analytic', 'color', 'red')
hold on
title('Combined intensity at $x=L=10^{-6}$m for excitation 1 and 2','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
legend()
drawnow


%% Calculate errors

% Calculate errors
error_x0 = norm((S_ex1(n_screen_0) + S_ex2(n_screen_0)) - (S1_f1 + S1_f2)) / norm(S1_f1 + S1_f2);
error_xL = norm((S_ex1(n_screen_L) + S_ex2(n_screen_L)) - (S3_f1 + S3_f2)) / norm(S3_f1 + S3_f2);

% Display errors
disp(['Relative error of intensity at screen at x = 0: ', num2str(error_x0)]);
disp(['Relative error of intensity at screen at x = L: ', num2str(error_xL)]);
