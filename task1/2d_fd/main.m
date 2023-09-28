%% Task1: Double slit

clc
clear all

% Path mesh functions
path_msh_func = '../fit/2d/mesh';
path_mat_func = '../fit/2d/matrices';
path_solver_func = '../fit/2d/solver';
path_solver_util = '../fit/2d/util';
path_util_func = '../fit/util';
path_verify_func = '../task1/verifications';
path_2d_fd_func = '../task1/2d_fd';

% Add paths
cd('../');
addpath(path_msh_func, path_mat_func, path_solver_func, path_solver_util, path_util_func, path_verify_func, path_2d_fd_func)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_farfield = 0;          % Calculate the fresnel number and test the farfield condition
use_y_symmetry = 1;         % Whether to use the symmetry in y direction
plot_mesh = 0;              % Plot the 2D mesh
solve_eq = 1;               % Solve the 2D Helmholtz equation
plot_field = 1;             % Plot the 2D electrical field
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
omega1 = 2*pi*f1;   % [1/s]
E1 = 250;           % [V/m]

lambda2 = 510e-9;   % [m]
f2 = c/lambda2;     % [Hz]
omega2 = 2*pi*f2;   % [1/s]
E2 = 500;           % [V/m]

if test_farfield
    fprintf('Fresnel number = %f for wave 1', fresnel_number(delta, L, lambda1))
    fprintf('Fresnel number = %f for wave 2', fresnel_number(delta, L, lambda2))
end

% Problem size in wavelength        |   
%               L2                  |      b: middle index
%                y                  |         
%       ###################         |      d: distance between the
%       #                 #         |         excitations in wave-
%  L3   ->                #   L1    |         length  
%  x    #      Model      #         |
%       ->                #         |      L: side index bc
%       #                 #         |         [L1, L2, L3, L4]
%       ###################         |
%              L4

d = 4e-6;                   % slit distance
delta = 1e-6;               % slit width
h = 8e-6;                   % screen height
L = 10e-6;                  % screen distance
NPML = [20, 20, 20, 20];    % [L1, L2, L3, L4]; 0,1:=PMC
offset = NPML;              % Total offset from boundaries

%% Generate Mesh
elem_per_wavelength = 20;

dx = lambda1*(offset(3) + offset(1))/elem_per_wavelength;  % Extra space in x direction for PML
xmesh = linspace(0, L + dx, ceil( (L + dx)/lambda1*elem_per_wavelength) );

if use_y_symmetry
    NPML(2) = 0;
    offset(2) = 0;
    dy = lambda1*(offset(4) + offset(2))/elem_per_wavelength;  % Extra space in y direction for PML
    ymesh = linspace(0, h/2 + dy, ceil( (h/2 + dy)/lambda1*elem_per_wavelength ));
else
    dy = lambda1*(offset(4) + offset(2))/elem_per_wavelength;  % Extra space in y direction for PML
    ymesh = linspace(-(h + dy)/2, (h + dy)/2, ceil( (h + dy)/lambda1*elem_per_wavelength ));
end

msh = cartMesh_2D(xmesh, ymesh);

% Set rhs and bc vectors
idx_bc = calc_slit_idx(msh, d, delta, use_y_symmetry) + NPML(3); % Transform y-indices to canonical index
jsbow = sparse(3*msh.np, 1);
ebow1_bc = NaN(3*msh.np, 1);
ebow2_bc = NaN(3*msh.np, 1);
ebow1_bc(idx_bc) = E1;
ebow2_bc(idx_bc) = E2;

if plot_mesh
    [X,Y] = meshgrid(xmesh, ymesh);
    plot(X, Y, 'blue', X', Y', 'blue')
    hold on
    plot(0, d/2, 'x', 'color', '#A2142F', 'linewidth', 6)
    plot(0, -d/2, 'x', 'color', '#A2142F', 'linewidth', 6)
    xline(0, 'red')
    xline(L, 'red')
    yline(-h/2, 'red')
    yline(h/2, 'red')
end


%% Solve Helmholtz

% Create matrices
[C, ~, ~] = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);

meps = createMeps_2D(msh, ds, da, dat, ones(msh.np, 1), eps);
mmui = createMmui_2D(msh, ds, dst, da, ones(msh.np, 1), mui);

% Solve
if solve_eq
    ebow1 = frequency_domain_2D(msh, C, meps, mmui, 0, jsbow, idx_bc, ebow1_bc(idx_bc), omega1, NPML);
    ebow2 = frequency_domain_2D(msh, C, meps, mmui, 0, jsbow, idx_bc, ebow2_bc(idx_bc), omega2, NPML);
    ebow = ebow1 + ebow2;
end


%% Postprocessing
idx = 1+offset(4):length(ymesh)-offset(2);  % Indices at which to evaluate the field
y = ymesh(idx);                             % y values at those indices

if plot_field
    figure
    [X,Y] = meshgrid(msh.xmesh, msh.ymesh);
    e_surf = reshape(real(ebow(2*msh.np+1:3*msh.np)), [msh.nx, msh.ny]);
    e_surf_plot = surf(X,Y,e_surf');
    %xlim([0, L])
    ylim([-h/2, h/2])
    set(e_surf_plot,'LineStyle','none')
    set(gca,'ColorScale','log')
end

% Intensity calculation % TODO: CAN'T add intensities!!!
idx_screen = msh.nx * (1:msh.ny) - offset(1) + 2*msh.np;
e1_screen = ebow1(idx_screen)';
e2_screen = ebow2(idx_screen)';
%e_screen = ebow(msh.nx * (1:msh.ny) - NPML(1))';
e1_screen = e1_screen(idx);
e2_screen = e2_screen(idx);
%e_screen = e_screen(1+NPML(2):end-NPML(4));
I1 = c*eps/2 * abs(e1_screen).^2;
I2 = c*eps/2 * abs(e2_screen).^2;
I = I1 + I2;

if plot_intensity
    figure
    plot(y, I/max(I), 'DisplayName', 'Numerical', 'color', '#1e8080')
    hold on
    title('Intensity at the screen at $x=L=10^6$m','Interpreter','latex')
    xlabel('Position at the screen $y$ (m)','Interpreter','latex')
    ylabel('Intensity $I$','Interpreter','latex')
    xlim([-h/2, h/2])
    legend()
end

if plot_intensity_colored
    plot(y, I1/max(I), 'DisplayName', 'Wave 1', 'color', '#3d00ff')
    plot(y, I2/max(I), 'DisplayName', 'Wave 2', 'color', '#00ff00')
end


%% verifications

% Double slit and helmholtz formula
I1_farfield = intensity_farfield(E1, lambda1, d, delta, L, y);
I2_farfield = intensity_farfield(E2, lambda2, d, delta, L, y);
I_farfield = I1_farfield + I2_farfield;
I1_helmholtz = intensity_helmholtz(E1, lambda1, d, delta, L, y, ceil(length(idx_bc)/2));
I2_helmholtz = intensity_helmholtz(E2, lambda2, d, delta, L, y, ceil(length(idx_bc)/2));
I_helmholtz = I1_helmholtz + I2_helmholtz;
if plot_intensity_ana % TODO: FIX normalization
    plot(y, I_farfield/max(I_farfield), 'r--', 'DisplayName', 'Analytical (farfield)')
    plot(y, I_helmholtz/max(I_helmholtz), 'b--', 'DisplayName', 'Analytical (Helmholtz)')
end

% Error calculation
I_err = norm(I/max(I) - I_farfield/max(I_farfield))/norm(I_farfield/max(I_farfield));
I_err_helmholtz = norm(I/max(I) - I_helmholtz/max(I_helmholtz))/norm(I_helmholtz/max(I_helmholtz));
if calc_intensity_err
    fprintf('Relative L2 error between numerical and farfield solution = %f%% \n', 100*I_err)
    fprintf('Relative L2 error between numerical and Helmholtz solution = %f%%', 100*I_err_helmholtz)
end
