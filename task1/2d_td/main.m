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
use_y_symmetry = 0;         % Whether to use the symmetry in y direction
polarisation = 'z';         % Direction of polarisation of the electric field
plot_field = 1;             % Plot the 2D electrical field
plot_intensity = 1;         % Plot the numerically calculated intensity on the screen
plot_intensity_colored = 0; % Plot the calculated intensities in the actual light colors
plot_intensity_ana = 1;     % Plot the analytically calculated intensity on the screen
calc_intensity_err = 1;     % Calculates the error between analytical and numerical solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Problem Definition
c = 3e8;            % [m/s]

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

offset = [20, 20, 20, 20];                  % Total offset from boundaries
bc.bc = ["OPEN", "OPEN", "OPEN", "OPEN"];   % [L1, L2, L3, L4]
bc.NPML = offset;

%% Generate Mesh

% TD params
dt = 8e-17;
tend = 20/f1;
nt = ceil(tend/dt);

% Geo params
elem_per_wavelength = 10;
dx = lambda1*(offset(2) + offset(4))/elem_per_wavelength;  % Extra space in x direction for PML
xmesh = linspace(0, L + dx, ceil( (L + dx)/lambda1*elem_per_wavelength) );

if use_y_symmetry
    offset(1) = 0;
    bc.NPML(1) = 0;
    bc.bc(1) = 'PMC';
    dy = lambda1*(offset(1) + offset(3))/elem_per_wavelength;  % Extra space in y direction for PML
    ymesh = linspace(0, h/2 + dy, ceil( (h/2 + dy)/lambda1*elem_per_wavelength ));
else
    dy = lambda1*(offset(1) + offset(3))/elem_per_wavelength;  % Extra space in y direction for PML
    ymesh = linspace(-(h + dy)/2, (h + dy)/2, ceil( (h + dy)/lambda1*elem_per_wavelength ));
end

msh = cartMesh_2D(xmesh, ymesh);

% Material params
material_regions.epsilon0 = 8.854187e-12;
material_regions.mu0i = 1/(pi*4e-7);

boxesEpsilonR(1).box = [1, msh.nx, 1, msh.ny];
boxesEpsilonR(1).value = 1;
material_regions.boxesEpsilonR = boxesEpsilonR;

boxesMuiR(1).box = [1, msh.nx, 1, msh.ny];
boxesMuiR(1).value = 1;
material_regions.boxesMuiR = boxesMuiR;

% Set rhs and bc vectors
idx_bc = calc_slit_idx(msh, d, delta, use_y_symmetry, polarisation) + offset(4); % Transform y-indices to canonical index
jsbow_bc = NaN(3*msh.np, 1);
ebow_bc = NaN(3*msh.np, 1);
ebow_bc(idx_bc) = 1;


%% Solve Helmholtz

% Apply bc and create matrices
[bc, W, ebow_bc, jsbow] = apply_bc(msh, bc, ebow_bc, jsbow_bc);
[MAT, bc] = generate_MAT(msh, bc, material_regions, ["CurlP"]);
MAT.mepsi = nullInv(MAT.meps);
[mur_edges,mur_n_edges, mur_deltas] = initMur_2D(msh, bc);

% Excitation
ebow_exi = @(t) ebow_bc * (E1*cos(2*pi*f1*t) + E2*cos(2*pi*f2*t));

% Vectors
ebow = sparse(3*msh.np, nt);
hbow = sparse(3*msh.np, nt);
ebow_abs = sparse(msh.np, nt);

% Solve
for i = 1:nt
    t = (i-1)*dt;

    % Save old and new values
    ebow_old = ebow(:,i);
    hbow_old = hbow(:,i);
    e_exi_old = ebow_exi(t);
    e_exi_new = ebow_exi(t+1);

    % Calc leapfrog
    [ebow_new, hbow_new] = solve_leapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.c,dt,W);

    % Apply open boundary with mur cond
    ebow_new = applyMur_2D(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt, bc);

    % Save ebow and hbow
    ebow(:,i+1) = ebow_new;
    hbow(:,i+1) = hbow_new;
    ebow_abs(:,i+1) = calc_abs_field(msh,ebow_new);

    if plot_field & mod(i, 10)
        [X,Y] = meshgrid(msh.xmesh, msh.ymesh);
        e_surf = reshape(ebow_abs(:,i+1), [msh.nx, msh.ny]);

        figure(1)
        e_surf_plot = surf(X,Y,e_surf');

        xlim([0, L])
        ylim([-h/2, h/2])
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

%% Postprocessing

return

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
    xlim([-h/2, h/2])
    legend()
end

if plot_intensity_colored
    plot(y, I1, 'DisplayName', 'Wave 1', 'color', '#3d00ff')
    plot(y, I2, 'DisplayName', 'Wave 2', 'color', '#00ff00')
end


%% verifications

% Double slit and helmholtz formula
I1_farfield = intensity_farfield(E1, lambda1, d, delta, L, y);
I2_farfield = intensity_farfield(E2, lambda2, d, delta, L, y);
I_farfield = I1_farfield + I2_farfield;
I_farfield = I_farfield/max(I_farfield);

I1_helmholtz = intensity_helmholtz(E1, lambda1, d, delta, L, y, ceil(length(idx_bc)/2));
I2_helmholtz = intensity_helmholtz(E2, lambda2, d, delta, L, y, ceil(length(idx_bc)/2));
I_helmholtz = I1_helmholtz + I2_helmholtz;
I_helmholtz = I_helmholtz/max(I_helmholtz);

if plot_intensity_ana
    plot(y, I_farfield, 'r--', 'DisplayName', 'Analytical (farfield)')
    plot(y, I_helmholtz, 'b--', 'DisplayName', 'Analytical (Helmholtz)')
end

% Error calculation
I_err = norm(I - I_farfield)/norm(I_farfield);
I_err_helmholtz = norm(I - I_helmholtz)/norm(I_helmholtz);
if calc_intensity_err
    fprintf('Relative L2 error between numerical and farfield solution = %f%% \n', 100*I_err)
    fprintf('Relative L2 error between numerical and Helmholtz solution = %f%%', 100*I_err_helmholtz)
end
