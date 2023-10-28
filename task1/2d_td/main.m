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
path_2d_td_func = '../task1/2d_td';

% Add paths
cd('../');
addpath(path_msh_func, path_mat_func, path_solver_func, path_solver_util, path_util_func, path_verify_func, path_2d_td_func)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_farfield = 1;          % Calculate the fresnel number and test the farfield condition
use_y_symmetry = 0;         % Whether to use the symmetry in y direction
polarisation = 'z';         % Direction of polarisation of the electric field
plot_field = 0;             % Plot the 2D electrical field
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

if test_farfield
    fprintf('Fresnel number = %f for wave 1', fresnel_number(delta, L, lambda1))
    fprintf('Fresnel number = %f for wave 2', fresnel_number(delta, L, lambda2))
end

bc.bc = ["OPEN", "OPEN", "OPEN", "OPEN"];   % [L1, L2, L3, L4]

%% Generate Mesh

% Geo params
elem_per_wavelength = 10;
wavelengths_pml = 6;
bc.NPML = [1,1,1,1]*wavelengths_pml*elem_per_wavelength;
offset =  [1,1,1,1]*(wavelengths_pml + 1)*elem_per_wavelength;
dx = lambda1*(offset(2) + offset(4))/elem_per_wavelength;  % Extra space in x direction for offset
xmesh = linspace(0, L + dx, ceil( (L + dx)/lambda1*elem_per_wavelength) );

if use_y_symmetry
    offset(1) = 0;
    bc.bc(1) = 'PMC';
    dy = lambda1*offset(3)/elem_per_wavelength;  % Extra space in y direction for offset
    ymesh = linspace(0, h/2 + dy, ceil( (h/2 + dy)/lambda1*elem_per_wavelength ));
else
    dy = lambda1*(offset(1) + offset(3))/elem_per_wavelength;  % Extra space in y direction for offset
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

%% Apply bc and create matrices
[bc, W, ebow_bc, jsbow] = apply_bc(msh, bc, ebow_bc, jsbow_bc);
[MAT, bc] = generate_MAT(msh, bc, material_regions, ["CurlP"]);
MAT = conductivePML_2D(bc, msh, MAT, f1);
MAT.mepsi = nullInv(MAT.meps);

% TD params
dt = CFL(msh, MAT);
t_spread = sqrt((1/2*h)^2 + L^2)/c;
tend = t_spread + 5/f1;
nt = ceil(tend/dt);

% Excitation
ebow_exi_f1 = @(t) ebow_bc * (E1*cos(2*pi*f1*t));
ebow_exi_f2 = @(t) ebow_bc * (E2*cos(2*pi*f2*t));


%% Solve for f1 
% init Vectors
ebow = zeros(3*msh.np, nt);
hbow = zeros(3*msh.np, nt);
ebow_abs_f1 = zeros(msh.np, nt);

% Solve with leapfrog
for i = 1:nt
    t = (i-1)*dt;

    % Save old and new values
    ebow_old = ebow(:,i);
    hbow_old = hbow(:,i);
    e_exi_old = ebow_exi_f1(t);
    e_exi_new = ebow_exi_f1(t+1);

    % Calc leapfrog
    [ebow_new,hbow_new] = solve_FullLeapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.kaps,MAT.c,dt,W);

    % Save ebow and hbow
    ebow(:,i+1) = ebow_new;
    hbow(:,i+1) = hbow_new;
    ebow_abs_f1(:,i+1) = calc_abs_field(msh,ebow_new);

    if plot_field & mod(i, 10)
        [X,Y] = meshgrid(msh.xmesh, msh.ymesh);
        e_surf = reshape(ebow_abs_f1(:,i+1), [msh.nx, msh.ny]);

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

%% Solve for f2
% init Vectors
ebow = zeros(3*msh.np, nt);
hbow = zeros(3*msh.np, nt);
ebow_abs_f2 = zeros(msh.np, nt);

% Solve with leapfrog
for i = 1:nt
    t = (i-1)*dt;

    % Save old and new values
    ebow_old = ebow(:,i);
    hbow_old = hbow(:,i);
    e_exi_old = ebow_exi_f2(t);
    e_exi_new = ebow_exi_f2(t+1);

    % Calc leapfrog
    [ebow_new,hbow_new] = solve_FullLeapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.kaps,MAT.c,dt,W);

    % Save ebow and hbow
    ebow(:,i+1) = ebow_new;
    hbow(:,i+1) = hbow_new;
    ebow_abs_f2(:,i+1) = calc_abs_field(msh,ebow_new);

    if plot_field & mod(i, 10)
        [X,Y] = meshgrid(msh.xmesh, msh.ymesh);
        e_surf = reshape(ebow_abs_f2(:,i+1), [msh.nx, msh.ny]);

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
% Calculate Intensities
[I_f1,y] = calc_intensity(msh, ebow_abs_f1, offset);
[I_f2,y] = calc_intensity(msh, ebow_abs_f2, offset);
I_sum = I_f1 + I_f2;

% Norm Intensities
if max(I_f1) ~= 0
    I_f1 = I_f1/max(I_f1);
end

if max(I_f2) ~= 0
    I_f2 = I_f2/max(I_f2);
end

if max(I_sum) ~= 0
    I_sum = I_sum/max(I_sum);
end

% Calculate min and max locations
[d_max_f1,d_min_f1] = calc_max_min_pos(y, L, d, lambda1);
[d_max_f2,d_min_f2] = calc_max_min_pos(y, L, d, lambda2);

%% componentwise analytically solution (farfield)
% calc farfield analytical solution
I1_farfield = intensity_farfield(E1, lambda1, d, delta, L, y);
I2_farfield = intensity_farfield(E2, lambda2, d, delta, L, y);
Isum_farfield = I1_farfield + I2_farfield;
% norm intensities
I1_farfield = I1_farfield / max(I1_farfield);
I2_farfield = I2_farfield / max(I2_farfield);
Isum_farfield = Isum_farfield / max(Isum_farfield);

figure
plot(y, I1_farfield, 'DisplayName', 'Farfield f1', 'color', '#77AC30')
hold on
plot(y, I2_farfield, 'DisplayName', 'Farfield f2', 'color', '#3d00ff')
hold on
scatter(d_max_f1, zeros(max(size(d_max_f1)),1),'filled', 'DisplayName', 'max f1', 'color', 'green');
hold on;
scatter(d_min_f1, zeros(max(size(d_min_f1)),1),'DisplayName', 'min f1', 'color', 'green');
hold on
scatter(d_max_f2, zeros(max(size(d_max_f2)),1),'filled', 'DisplayName', 'max f2', 'color', 'blue');
hold on;
scatter(d_min_f2, zeros(max(size(d_min_f2)),1),'DisplayName', 'min f2', 'color', 'blue');
title('Analytical solution (farfield) at the screen at $x=L=10^6$m','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
xlim([-h/2, h/2])
legend()

%% componentwise analytically solution (helmholtz)
% calc helmholtz analytical solution
I1_helmholtz = intensity_helmholtz(E1, lambda1, d, delta, L, y, ceil(length(idx_bc)/2));
I2_helmholtz = intensity_helmholtz(E2, lambda2, d, delta, L, y, ceil(length(idx_bc)/2));
Isum_helmholtz = I1_helmholtz + I2_helmholtz;
% norm intensities
I1_helmholtz = I1_helmholtz / max(I1_helmholtz);
I2_helmholtz = I2_helmholtz / max(I2_helmholtz);
Isum_helmholtz = Isum_helmholtz / max(Isum_helmholtz);

figure
plot(y, I1_helmholtz, 'DisplayName', 'Helmholtz f1', 'color', '#77AC30')
hold on
plot(y, I2_helmholtz, 'DisplayName', 'Helmholtz f2', 'color', '#3d00ff')
title('Analytical solution (helmoltz) at the screen at $x=L=10^6$m','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
xlim([-h/2, h/2])
legend()

%% componentwise comparison (for f1)
figure
plot(y, I1_farfield, 'DisplayName', 'Farfield f1', 'color', '#77AC30')
hold on
plot(y, I1_helmholtz, 'DisplayName', 'Helmholtz f1', 'color', '#3d00ff')
hold on
plot(y, I_f1, 'DisplayName', 'Numeric f1', 'color', '#D95319')
title('Comparison Farfield, Helmholtz and Numeric solution for f1','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
xlim([-h/2, h/2])
legend()

%% combined wave comparison
figure
plot(y, I_sum, 'DisplayName', 'Numerical', 'color', '#D95319')
hold on
plot(y, Isum_helmholtz, 'DisplayName', 'Helmholtz', 'color', '#3d00ff')
hold on
plot(y, Isum_farfield, 'DisplayName', 'Farfield', 'color', '#77AC30')
title('Intensity combined waves at $x=L=10^6$m','Interpreter','latex')
xlabel('Position at the screen $y$ (m)','Interpreter','latex')
ylabel('Intensity $I$','Interpreter','latex')
xlim([-h/2, h/2])
legend()

%% Error calculation
I_err_f1    = norm(I_f1 - I1_farfield)/norm(I1_farfield);
I_err_f2    = norm(I_f2 - I2_farfield)/norm(I2_farfield);
I_err       = norm(I_sum - Isum_farfield)/norm(Isum_farfield);

I_err_f1_helmholtz  = norm(I_f1 - I1_helmholtz)/norm(I1_helmholtz);
I_err_f2_helmholtz  = norm(I_f2 - I2_helmholtz)/norm(I2_helmholtz);
I_err_helmholtz     = norm(I_sum - Isum_helmholtz)/norm(Isum_helmholtz);
if calc_intensity_err
    fprintf('Rel. L2 err numerical and farfield f1 = %f%% \n', 100*I_err_f1)
    fprintf('Rel. L2 err numerical and farfield f2 = %f%% \n', 100*I_err_f2)
    fprintf('Rel. L2 err numerical and farfield combined = %f%% \n', 100*I_err)

    fprintf('Rel. L2 err numerical and helmholtz f1 = %f%% \n', 100*I_err_f1_helmholtz)
    fprintf('Rel. L2 err numerical and helmholtz f2 = %f%% \n', 100*I_err_f2_helmholtz)
    fprintf('Rel. L2 err numerical and helmholtz combined = %f%% \n', 100*I_err_helmholtz)
end
