%% Task2: Slit with reflection

clc
clear all

% Path mesh functions
path_msh_func = '../fit/2d/mesh';
path_mat_func = '../fit/2d/matrices';
path_solver_func = '../fit/2d/solver';
path_solver_util = '../fit/2d/util';
path_util_func = '../fit/util';
path_verify_func = '../task2/verifications';
path_2d_fd_func = '../task2/2d_fd';

% Add paths
cd('../');
addpath(path_msh_func, path_mat_func, path_solver_func, path_solver_util, path_util_func, path_verify_func, path_2d_fd_func)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_farfield = 0;          % Calculate the fresnel number and test the farfield condition
polarisation = 'z';         % Direction of polarisation of the electric field
solve_eq = 1;               % Solve the 2D Helmholtz equation
plot_field = 0;             % Plot the 2D electrical field
plot_intensity = 1;         % Plot the numerically calculated intensity on the screen
plot_intensity_colored = 0; % Plot the calculated intensities in the actual light colors
plot_intensity_ana = 1;     % Plot the analytically calculated intensity on the screen
calc_intensity_err = 1;     % Calculates the error between analytical and numerical solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Problem Definition
c = 3e8;            % [m/s]
eps = 8.854e-12;
mui = 1/(4*pi*1e-7);

lambda1 = 430e-9;   % [m]
f1 = c/lambda1;     % [Hz]
E1 = 250;           % [V/m]

lambda2 = 510e-9;   % [m]
f2 = c/lambda2;     % [Hz]
E2 = 500;           % [V/m]

if test_farfield
    fprintf('Fresnel number = %f for wave 1', fresnel_number(delta, L, lambda1))
    fprintf('Fresnel number = %f for wave 2', fresnel_number(delta, L, lambda2))
end

% Problem size in wavelength        |
%               L3                  |      b: middle index
%                y                  |
%       ###################         |      d: distance between the
%       #                 #         |         excitations in wave-
%  L4   ->                #   L2    |         length
%  x    #      Model      #         |
%       ->                #         |      L: side index bc
%       #                 #         |         [L1, L2, L3, L4]
%       ###################         |
%              L1

d = 4e-6;       % slit distance
delta = 1e-6;   % slit width
h = 8e-6;       % screen height
L = 10e-6;      % screen distance

offset = [0, 20, 20, 20];                   % Total offset from boundaries
bc.bc = ["PEC", "OPEN", "OPEN", "OPEN"];    % [L1, L2, L3, L4]
bc.NPML = offset;

%% Generate Mesh
elem_per_wavelength = 15;
dx = lambda1*(offset(2) + offset(4))/elem_per_wavelength;  % Extra space in x direction for PML
xmesh = linspace(0, L + dx, ceil( (L + dx)/lambda1*elem_per_wavelength) );
dy = lambda1*(offset(1) + offset(3))/elem_per_wavelength;  % Extra space in y direction for PML
ymesh = linspace(0, h/2 + dy, ceil( (h/2 + dy)/lambda1*elem_per_wavelength ));
msh = cartMesh_2D(xmesh, ymesh);

% Set rhs and bc vectors
idx_bc = calc_slit_idx(msh, d, delta, polarisation) + offset(4); % Transform y-indices to canonical index
jsbow = sparse(3*msh.np, 1);
ebow1_bc = NaN(3*msh.np, 1);
ebow2_bc = NaN(3*msh.np, 1);
ebow1_bc(idx_bc) = E1;
ebow2_bc(idx_bc) = E2;


%% Solve Helmholtz

% Create matrices
C = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);

meps = createMeps_2D(msh, ds, da, dat, ones(msh.np, 1), eps);
mmui = createMmui_2D(msh, ds, dst, da, ones(msh.np, 1), mui);

% Solve
if solve_eq
    [bc, W, e_exi, jsbow] = apply_bc(msh, bc, ebow1_bc, jsbow);
    ebow1 = solve_helmholtz_2d_fd(msh, W, C, meps, mmui, jsbow, e_exi, f1, bc);
    [bc, W, e_exi, jsbow] = apply_bc(msh, bc, ebow2_bc, jsbow);
    ebow2 = solve_helmholtz_2d_fd(msh, W, C, meps, mmui, jsbow, e_exi, f2, bc);
    ebow = ebow1 + ebow2;
end


%% Postprocessing

% Field plot
ebow_abs = calc_abs_field(msh,ebow);

if plot_field
    figure
    [X,Y] = meshgrid(msh.xmesh, msh.ymesh);
    e_surf = reshape(ebow_abs, [msh.nx, msh.ny]);
    e_surf_plot = surf(X,Y,e_surf');
    xlim([0, L])
    ylim([0, h/2])
    set(e_surf_plot,'LineStyle','none')
    colormap hot;
    title('Absolute value of electric field','Interpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex')
    ylabel('$y$ (m)','Interpreter','latex')
    zlabel('Absolute value','Interpreter','latex')
end

% Intensity calculation
[I1,y] = calc_intensity(msh, ebow1', offset);
I2 = calc_intensity(msh, ebow2', offset);
I = I1 + I2;
I1 = I1/max(I);
I2 = I2/max(I);
I = I/max(I);

if plot_intensity
    figure
    plot(y, I, 'DisplayName', 'Numerical', 'color', '#1e8080')
    hold on
    title('Intensity at the screen at $x=L=10^6$m','Interpreter','latex')
    xlabel('Position at the screen $y$ (m)','Interpreter','latex')
    ylabel('Intensity $I$','Interpreter','latex')
    xlim([0, h/2])
    legend()
end

if plot_intensity_colored
    plot(y, I1, 'DisplayName', 'Wave 1', 'color', '#3d00ff')
    plot(y, I2, 'DisplayName', 'Wave 2', 'color', '#00ff00')
end


%% verifications

% Analytical helmholtz formula
I1_helmholtz = intensity_helmholtz(E1, lambda1, d, delta, L, y, ceil(length(idx_bc)/2));
I2_helmholtz = intensity_helmholtz(E2, lambda2, d, delta, L, y, ceil(length(idx_bc)/2));
I_helmholtz = I1_helmholtz + I2_helmholtz;
I_helmholtz = I_helmholtz/max(I_helmholtz);

if plot_intensity_ana
    plot(y, I_helmholtz, 'b--', 'DisplayName', 'Analytical (Helmholtz)')
end

% Error calculation
I_err = norm(I - I_helmholtz)/norm(I_helmholtz);
if calc_intensity_err
    fprintf('Relative L2 error between numerical and Helmholtz solution = %f%%', 100*I_err)
end
