%% Task3 - Thin film - time domain solver

% Clear variables
clc
clear

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
path_verifications = append(parent, '\task3\verifications');

% Add paths
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func, path_verifications)

%% Define excitation

% Add basic constants to material_regions object
material_regions.epsilon0 = 8.854187e-12;
material_regions.mu0i = 1/(pi*4e-7);

% Excitation 1
lambda_1    = 430e-9;
f1          = sqrt(material_regions.mu0i/material_regions.epsilon0)/lambda_1;           
E1          = 250;
% Excitation 2
lambda_2    = 510e-9;
f2          = sqrt(material_regions.mu0i/material_regions.epsilon0)/lambda_2;           
E2          = 500;
func_exi_1  = @(t)(E1 * sin(2*pi*f1*t)+ E2 * sin(2*pi*f2*t));


%% Define important parameters for the simulation

% intnsity plots ?
plot_intensity = 1;

% plot analytic ?
plot_analytic_ypol = 1;

% Polarization: 1 for z and 2 for y
polarization = 2;

% Elements per wavelength
elem_per_wavelength = 10;

% Offset in each direction
offset = [0,2,0,2];

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
x_basic = linspace(0, L * 1e-6, points_x);
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


%% Find important indices in xmesh and ymesh

% Index of x = 0
idx_x0      = find(xmesh == 0);

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
    n_exi = 1 + (idx_x0 - 1) * Mx + ((idx_ymh_h:idx_yph_h) - 1) * My;

    % Use corresponding edges in z-direction for excitation
    e_exitation(n_exi + 2*np) = 1;

% Generate excitation for polarization in y-direction
elseif polarization == 2

    % Determine indices for points in the single slit
    n_exi = 1 + (idx_x0 - 1) * Mx + ((idx_ymh_h:idx_yph_h-1) - 1) * My;

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


%% Set up parameters for the simulation in time domain for excitation 1 

% Time step size
dt = 4e-17;

% Fit dt to period of excitation 1
dt = 1/f1 / ceil(1/f1 / dt);

% Calculate end time regarding the time needed for the calculation of the
% on avarage emitted power
t_end = L * 1e-6 / sqrt(material_regions.mu0i/material_regions.epsilon0);
t_end = t_end + 2 * max(1/f1, 1/f2);


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

% Calculate time steps
steps_movie = 0;
for t = linspace(0,t_end,ceil(t_end/dt))

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
        S_ex1 = S_ex1 + CalcPowerSurfaceXY(msh, ebow_new, hbow_new, MAT.ds, MAT.dst, MAT.da);
        i_steps = i_steps + 1;
    end

    % Plot excitation of the electric field
    steps_movie = steps_movie + 1;

    if steps_movie == 10

        steps_movie = 0;

        [X,Y] = meshgrid(xmesh, ymesh);
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
        drawnow
        
    end

end

% Get time avarage of the power through each surface
S_ex1 = S_ex1/i_steps;


%% Plot avarage power on x=0 and x=L screen for excitation 1

if plot_intensity

    % solve analytic
    if polarization == 2
        % y-polarization solution at x=0 and x=L
        x_pos = [idx_x0, idx_xL]; 
        [~,~,analytic] = analytic_sol_ypol(x_pos, E1, E2, lambda_1, lambda_2, L*1e-6/2, a*1e-6, 2, le^2);
    else
        % ToDo z-polarization analytic = 0 (?)
        analytic = [0, 0];
    end

    % Calculate indices of screen at x = 0
    n_screen_0 = 1 + (idx_x0 - 1) * Mx + ((idx_ymh_h:idx_yph_h)-1) * My;
    % Calculate y-coordinates of corresponding faces as parts of the screen
    y_coord_screen = msh.ymesh + le/2;
    % Plot avarage power on screen over associated y-coordinates
    figure(2)
    plot(y_coord_screen, S_ex1(n_screen_0 + msh.np),'DisplayName', 'Numerical', 'color', '#1e8080')
    hold on
    plot(y_coord_screen, ones(1, length(y_coord_screen))*analytic(1), 'DisplayName', 'analytic', 'color', 'red')
    hold on
    title('Intensity at the screen at $x=0$m','Interpreter','latex')
    xlabel('Position at the screen $y$ (m)','Interpreter','latex')
    ylabel('Intensity $I$','Interpreter','latex')
    xlim([-h/2*1e-6, h/2*1e-6])
    ylim([0, 2e-5])
    legend()

    % Calculate indices of the second screen at x = L
    n_screen_L = 1 + (idx_xL - 1) * Mx + ((idx_ymh_h:idx_yph_h)-1) * My;
    % Plot avarage power on screen x=L over associated y-coordinates
    figure(3)
    plot(y_coord_screen, S_ex1(n_screen_L + msh.np), 'DisplayName', 'Numerical', 'color', '#1e8080')
    hold on
    plot(y_coord_screen, ones(1,length(y_coord_screen))*analytic(2), 'DisplayName', 'analytic', 'color', 'red')
    hold on
    title('Intensity at the screen at $x=L=10^6$m','Interpreter','latex')
    xlabel('Position at the screen $y$ (m)','Interpreter','latex')
    ylabel('Intensity $I$','Interpreter','latex')
    xlim([-h/2*1e-6, h/2*1e-6])
    ylim([0, 5e-5])
    legend()
end

if plot_analytic_ypol
    % Calculate analytic intensity across whole domain
    x_pos = msh.xmesh;
    [e_pos, h_pos, s_pos] = analytic_sol_ypol(x_pos, E1, E2, lambda_1, lambda_2, L*1e-6/2, a*1e-6, 2, le^2);
    % Plot avarage power on screen x=L over associated y-coordinates
    figure(4)
    plot(x_pos, e_pos./max(e_pos), 'DisplayName', 'electric field', 'color', '#1e8080')
    hold on
    plot(x_pos, s_pos./max(s_pos), 'DisplayName', 'intensity', 'color', 'red')
    hold on
    title('normalized analytic solutions across whole domain','Interpreter','latex')
    xlabel('Position in the domain $x$ (m)','Interpreter','latex')
    ylabel('Intensity $I$','Interpreter','latex')
    legend()
end

figure(5)
s_surf = reshape(S_ex1(idx_2_plot,1), [msh.nx, msh.ny]);
s_surf_plot = surf(X,Y,s_surf');
title('numeric calculated intensity','Interpreter','latex')
%xlim([0, L*1e-6])
%ylim([0, h/2*1e-6])
set(s_surf_plot,'LineStyle','none')
colormap hot;

