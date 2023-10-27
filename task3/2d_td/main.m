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
path_msh_func = '../fit/2d/mesh';
path_mat_func = '../fit/2d/matrices';
path_solver_func = '../fit/2d/solver';
path_util_func = '../fit/2d/util';

% Add paths
cd('../');
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

func_exi =@(t) func_exi_1(t) + func_exi_2(t);

%% Define important parameters for the simulation

% Polarization: 1 for z and 2 for y
polarization = 2;

% show field plots during calculation
field_plots = true;

% mesh refinement in the thin film
mesh_refinement = false;

% Elements per wavelength
elem_per_wavelength = 18;

% Offset in each direction
offset = [0,3*elem_per_wavelength,0,6*elem_per_wavelength];

% Edit boundary conditions
if polarization == 1
    bc.bc = ["PMC", "OPEN", "PMC", "OPEN"];
elseif polarization == 2
    bc.bc = ["PEC", "OPEN", "PEC", "OPEN"];
end


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

% % Calculate offset in micro meter
% offset = offset * num_e;

% Calculate number of points in x- and y-direction
points_x = num_e * L + 1;
points_y = num_e * h + 1;

% Calculate xmesh with respect to the choosen offset
x_offset1 = (-offset(4):-1) * le;
x_offset2 = L * 1e-6 + (1:offset(2)) * le;
% refinement of factor 4 in thin film
if mesh_refinement == true
    x_basic = [linspace(0, (L/2)*1e-6, num_e*(L/2)+1), ...
        linspace((L/2)*1e-6 + le/2, (L/2+a)*1e-6 - le/2, 2*num_e*a-2), ...
        linspace((L/2+a)*1e-6, L*1e-6, num_e*(L/2-a)+1)];
else
    x_basic = [linspace(0, (L/2) * 1e-6, num_e*(L/2)+1), linspace((L/2) * 1e-6 + le, (L) * 1e-6, num_e*(L/2))];
end
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
idx_x0      = find(xmesh == 0);

% Index in x-direction for excitation (inside the offset!!!!)
idx_exi     = idx_x0-(offset(4)/2);

% Index of x = L/2
idx_xL_h    = find(xmesh == L/2 * 1e-6);

% Index of x = L/2 + a
idx_xL_h_a  = find(xmesh == (L/2 + a)* 1e-6);

% Index of x = L
idx_xL      = find(xmesh == L * 1e-6);

% Index of y = -h/2
idx_ymh_h   = find(ymesh == -h/2 * 1e-6);

% Index of y = 0
idx_y0      = find(ymesh == 0);

% Index of y = h/2
idx_yph_h   = find(ymesh == h/2 * 1e-6);


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
% Relative permittivity everywhere equal to one
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

% calc reflection arrival at x=0
c0 = sqrt(MAT.mu0i/MAT.epsilon0);
% time of wave arrival + four percent (for possible num. errors)
t_end = (offset(4)*le/2 + L*1e-6)/c0 *1.04


%% Set up parameters for the simulation in time domain for excitation 1 

% Time step size
dt = CFL(msh, MAT)* 0.7;

% Fit dt to period of excitation 1
%dt = 1/f1 / ceil(1/f2 / dt);

% Calculate end time regarding the time needed for the calculation of the
% on avarage emitted power
%t_end = (30+0.7)/f1;


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

% Indices for plot of electric field
if polarization == 1
    idx_2_plot = 2*np+1:3*np;
elseif polarization == 2
    idx_2_plot = np+1:2*np;
end

% init gridpoints for the plots
[X,Y] = meshgrid(msh.xmesh, msh.ymesh);

% Calculate time steps
steps_movie = 0;
for t = linspace(0,t_end,ceil(t_end/dt))

    % Old values
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    % Calculate value for excitation
    e_exi_old = e_exi * func_exi(t);
    e_exi_new = e_exi * func_exi(t+1);

    % Execute timestep with leapfrog
    [ebow_new,hbow_new] = solve_leapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.c,dt,W);

    % Apply open boundary with mur condition
    ebow_new = applyMur_2D(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt, bc);

    % Sum up power through each surface for one period before the end
    if t >= t_end - 1/f1
        S_ex1 = S_ex1 + CalcPoyntinvectorXY(msh, ebow_new, hbow_new, MAT.ds, MAT.dst);
        i_steps = i_steps + 1;
    end

    % Plot excitation of the electric field
    steps_movie = steps_movie + 1;

    if steps_movie == 10 && field_plots

        steps_movie = 0;
        field_surf_plot(msh,X,Y,idx_2_plot, ebow_new, L, h)
        drawnow
    end
end

if ~(field_plots)
    % plot only final ebow
    field_surf_plot(msh,X,Y,idx_2_plot, ebow_new, L, h)
end

% Get time avarage of the power through each surface
S_ex1 = S_ex1/i_steps;


%% Calculate analytical solution for the intensity for excitation 1
[S1_f1, S3_f1] = AnaSolPoyntin(E1, f1, 1 * material_regions.epsilon0, 4 * material_regions.epsilon0, material_regions.mu0i, material_regions.mu0i, a*1e-6);
[S1_f2, S3_f2] = AnaSolPoyntin(E2, f2, 1 * material_regions.epsilon0, 4 * material_regions.epsilon0, material_regions.mu0i, material_regions.mu0i, a*1e-6);

% add intensities
S1 = S1_f1 + S1_f2;
S3 = S3_f1 + S3_f2;

%% Plot avarage power on x=0 and x=L screen
%if polarization == 2
% Calculate indices of screen at x = 0 (y-surfaces)
n_screen_0 = 1 + (idx_x0 - 1) * Mx + ((idx_ymh_h:idx_yph_h-1)-1) * My + msh.np;
% Calculate y-coordinates of corresponding faces as parts of the screen
y_coord_screen = msh.ymesh(idx_ymh_h:idx_yph_h-1) + le/2;
% calculate relative error
error_x0 = abs(mean(S_ex1(n_screen_0))-real(S1))/real(S1)*100;
% Plot avarage power on screen over associated y-coordinates
figure(2)
plot(y_coord_screen, S_ex1(n_screen_0),'DisplayName', 'Numerical', 'color', '#1e8080')
hold on
plot(y_coord_screen, ones(1, length(y_coord_screen))*real(S1), 'DisplayName', 'Analytic', 'color', 'red')
hold on
title('Intensity at $x=0$m','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
xlim([-h/2*1e-6, h/2*1e-6])
%ylim([0, 2e-5])
legend()
disp(['relative intensity error at x=0: ',num2str(error_x0),' %'])


% Calculate indices of the second screen at x = L (y-surfaces)
n_screen_L = 1 + (idx_xL - 1) * Mx + ((idx_ymh_h:idx_yph_h-1)-1) * My +msh.np;
% calculate relative error
error_xL = abs(mean(S_ex1(n_screen_L))-real(S3))/real(S3)*100;
% Plot avarage power on screen x=L over associated y-coordinates
figure(3)
plot(y_coord_screen, S_ex1(n_screen_L), 'DisplayName', 'Numerical', 'color', '#1e8080')
hold on
plot(y_coord_screen, ones(1,length(y_coord_screen))*real(S3), 'DisplayName', 'Analytic', 'color', 'red')
hold on
title('Intensity at $x=L=10^6$m','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
xlim([-h/2*1e-6, h/2*1e-6])
%ylim([0, 5e-5])
legend()
disp(['relative intensity error at x=L: ',num2str(error_xL),' %'])
%end 

% numeric intensity in y-surfaces (= wave intensity)
figure(5)
s_surf = reshape(S_ex1((1+np:2*np),1), [msh.nx, msh.ny]);
s_surf_plot = surf(X,Y,real(s_surf'));
title('Numeric intensity across domain','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
xlim([0, L*1e-6])
ylim([0, h/2*1e-6])
set(s_surf_plot,'LineStyle','none')
colormap winter;

function field_surf_plot(msh, X, Y, idx_2_plot, ebow_new, L, h)
        e_surf = reshape(ebow_new(idx_2_plot,1), [msh.nx, msh.ny]);
        
        figure(1)
        e_surf_plot = surf(X,Y,e_surf');
    
        xlim([0, L*1e-6])
        ylim([0, h/2*1e-6])
        set(e_surf_plot,'LineStyle','none')
        view(2)
        colormap hot;
        title('Absolute value of electric field','Interpreter','latex')
        xlabel('$x$ (m)','Interpreter','latex')
        ylabel('$y$ (m)','Interpreter','latex')
        zlabel('Absolute value','Interpreter','latex')
end
