%% Task1: Double slit

clc
clear all

% Path mesh functions
path_msh_func = './fit/mesh';
path_mat_func = './fit/matrices';
path_solver_func = './fit/solver';
path_util_func = './fit/util';
path_verify_func = './task1/verifications';

% Add paths
cd('../');
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func, path_verify_func)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_farfield = 0;          % Calculate the fresnel number and test the farfield condition
plot_mesh = 0;              % Plot the 2D mesh
solve_eq = 1;               % Solve the 2D Helmholtz equation
plot_field = 0;             % Plot the 2D electrical field
plot_intensity = 1;         % Plot the numerically calculated intensity on the screen
plot_intensity_colored = 1; % Plot the calculated intensity in the actual light colors
plot_intensity_ana = 0;     % Plot the analytically calculated intensity on the screen
plot_intensity_err = 0;     % Plot the error between analytical and numerical solutions
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

d = 4e-6;       % slit distance
delta = 1e-6;   % slit width
h = 8e-6;       % screen height
L = 10e-6;      % screen distance
NPML = [20, 20, 20, 20];  % [L1, L2, L3, L4]; 0,1:=PMC

%% Generate Mesh
elem_per_wavelength = 15;
dx = lambda1*(NPML(3) + NPML(1))/elem_per_wavelength;  % Extra space in x direction
dy = lambda1*(NPML(4) + NPML(2))/elem_per_wavelength;  % Extra space in y direction
xmesh = linspace(0, L + dx, ceil( (L + dx)/lambda1*elem_per_wavelength) );
ymesh = linspace(-(h + dy)/2, (h + dy)/2, ceil( (h + dy)/lambda1*elem_per_wavelength ));
msh = cartMesh2D(xmesh, ymesh);

% Calculate BC indices
y_slit = [(-d - delta)/2, (-d + delta)/2, (d - delta)/2, (d + delta)/2]; % y values of upper and lower slit.
for i = 1:length(y_slit)
    % Find y-index closest to actual y_slit value
    [~,y_idx(i)] = min(abs(msh.ymesh - y_slit(i)));
end
y_idx = [y_idx(1):y_idx(2), y_idx(3):y_idx(4)]; % Find all y-indices between slits

% Set rhs and bc vectors
idx = msh.nx * (y_idx-1) + NPML(3); % Transform y-indices to canonical index
jsbow = sparse(msh.np, 1);
ebow1_bc = NaN(msh.np, 1);
ebow2_bc = NaN(msh.np, 1);
ebow1_bc(idx) = E1;
ebow2_bc(idx) = E2;

if test_farfield
    fprintf('Fresnel number = %f for wave 1', fresnel_number(delta, L, lambda1))
    fprintf('Fresnel number = %f for wave 2', fresnel_number(delta, L, lambda2))
end

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


%% Solution
if solve_eq
    ebow1 = solveHelmholtzTE(msh, eps, mui, jsbow, ebow1_bc, omega1, NPML);
    ebow2 = solveHelmholtzTE(msh, eps, mui, jsbow, ebow2_bc, omega2, NPML);
    ebow = ebow1 + ebow2;
    %save('ebow.mat', 'ebow')
end


%% Postprocessing
if plot_field
    figure
    [X,Y] = meshgrid(msh.xmesh, msh.ymesh);
    e_surf = reshape(real(ebow), [msh.nx, msh.ny]);
    e_surf_plot = surf(X,Y,e_surf');
    set(e_surf_plot,'LineStyle','none')
    set(gca,'ColorScale','log')
end

% Intensity calculation
e1_screen = ebow1(msh.nx * (1:msh.ny) - NPML(1))';
e2_screen = ebow2(msh.nx * (1:msh.ny) - NPML(1))';
%e_screen = ebow(msh.nx * (1:msh.ny) - NPML(1))';
e1_screen = e1_screen(1+NPML(2):end-NPML(4));
e2_screen = e2_screen(1+NPML(2):end-NPML(4));
%e_screen = e_screen(1+NPML(2):end-NPML(4));
y = ymesh(1+NPML(2):end-NPML(4));
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
I1_ana = intensity_ana(E1, lambda1, d, delta, L, y);
I2_ana = intensity_ana(E2, lambda2, d, delta, L, y);
I_ana = I1_ana + I2_ana;
I1_ana_helmholtz = helmholtz_ana(E1, lambda1, d, delta, L, y, ceil(length(y_idx)/2));
I2_ana_helmholtz = helmholtz_ana(E2, lambda2, d, delta, L, y, ceil(length(y_idx)/2));
I_ana_helmholtz = I1_ana_helmholtz + I2_ana_helmholtz;
if plot_intensity_ana % TODO: FIX normalization
    plot(y, I_ana/max(I_ana), 'r--', 'DisplayName', 'Analytical (double slit)')
    plot(y, I_ana_helmholtz/max(I_ana_helmholtz), 'b--', 'DisplayName', 'Analytical (helmholtz)')
end

% Error calculation
I_err = abs(I - I_ana);
if plot_intensity_err
    plot(y, I_err, 'DisplayName', 'Absolute error')
end
