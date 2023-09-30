%% Task4: Diffraction

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_field = 1;
plot_intensity = 1;
plot_intensity_ana = 1;
calc_intensity_err = 1;


%% Problem Definition
c0 = 3e8;            % [m/s]
eps0 = 8.854e-12;
mu0i = 1/(4*pi*1e-7);
n = 2;

v0 = c0;             % [m/s]
v1 = c0/n;           % [m/s] given c = 1/sqrt(mu0 * eps0)
% refraction idx n=2 by multiplying mu0 or eps0 for thin film area by 4
%eps1 = 4*eps;
rel_mui = 1;


lambda1 = 430e-9;   % [m]
f1 = c0/lambda1;     % [Hz]
omega1 = 2*pi*f1;   % [1/s]
E1 = 250;           % [V/m]

lambda2 = 510e-9;   % [m]
f2 = c0/lambda2;     % [Hz]
omega2 = 2*pi*f2;   % [1/s]
E2 = 500;           % [V/m]

% Problem size in wavelength        |   
%               L3                  |      b: middle index
%              y=-h/2               |         
%       ###################         |      d: distance between the
%       #                 #         |         excitations in wave-
%  L4   #                 #   L2    |         length  
%  x=0  ->     Model      #         |
%       #                 #         |      L: side index bc
%       #                 #         |         [L1, L2, L3, L4]
%       ###################         |
%              L1

%% geometry
h = 40e-6;     % [m]
L = 10e-6;    % [m]
%a = 100e-9;   % [m]
delta = 1e-6; % [m]

NPML = [20, 20, 20, 20];  % [L1, L2, L3, L4]; 0,1:=PMC
bc.bc = ["OPEN", "OPEN", "OPEN", "OPEN"];   % Total offset from boundaries
bc.NPML = NPML;

%% mesh
elem_per_wavelength = 10;

dx = lambda1*(NPML(3) + NPML(1))/elem_per_wavelength;  % Extra space in x direction for PML
xmesh = linspace(0, L + dx, ceil( (L + dx)/lambda1*elem_per_wavelength) );
dy = lambda1*(NPML(4) + NPML(2))/elem_per_wavelength;  % Extra space in y direction for PML
ymesh = linspace(-(h + dy)/2, (h + dy)/2, ceil( (h + dy)/lambda1*elem_per_wavelength ));
msh = cartMesh_2D(xmesh, ymesh);

% border x indices of different permittivity == thin film borders
%border1_x_idx = round((L/2)/(L/msh.nx)); 
%border2_x_idx = border1_x_idx + round(a/(L/msh.nx));
%actual_thickness = (border2_x_idx-border1_x_idx)*(L/msh.nx)

%% excitation

y_slit = [(- delta)/2, (delta)/2];
for i = 1:length(y_slit)
    % Find y-index closest to actual y_slit value
    [~,y_idx(i)] = min(abs(msh.ymesh - y_slit(i)));
end
idx = msh.nx * (y_idx-1) + 2*msh.np;

use_y_symmetry = 1;
polarisation = 'z'; 

%indices = 1:msh.np;
%idx = 2*msh.np+indices(~(mod(indices,msh.nx)-31));  % idx of all L2 boundary elements
idx_bc = calc_slit_idx(msh, delta, use_y_symmetry, polarisation) + NPML(3); % Transform y-indices to canonical index
jsbow = sparse(3*msh.np, 1);
ebow1_bc = NaN(3*msh.np, 1);
ebow2_bc = NaN(3*msh.np, 1);
ebow1_bc(idx_bc) = E1;
ebow2_bc(idx_bc) = E2;

% Anregung erstmal in z-Richtung (wie Task1). entspricht Aufgabe d)

%% solve system
%ebow1 = solve_helmholtz_2d_te_fd(msh, eps, mui, jsbow, ebow1_bc, omega1, NPML);
%ebow2 = solve_helmholtz_2d_te_fd(msh, eps, mui, jsbow, ebow2_bc, omega2, NPML);

[c, s, st] = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);
%rel_eps = sparse(1,length(ds)/3)+4;
%rel_mui = sparse(1,length(ds)/3)+1;
meps = createMeps_2D(msh, ds, da, dat, ones(msh.np, 1), eps0);
mmui = createMmui_2D(msh, ds, dst, da, ones(msh.np, 1), mu0i);
%mkaps = sparse(size(ds,1));

%idx_dof = isnan(ebow1_bc);
%idx_bc = ~idx_dof;

%[ebow1, hbow1] = frequency_domain_2D(msh, c, meps, mmui, 0, jsbow, idx_bc, ebow1_bc(idx_bc), omega1, NPML);
%[ebow2, hbow2] = frequency_domain_2D(msh, c, meps, mmui, 0, jsbow, idx_bc, ebow2_bc(idx_bc), omega2, NPML);

[bc, W, e_exi, jsbow] = apply_bc(msh, bc, ebow1_bc, jsbow);
ebow1 = solve_helmholtz_2d_fd(msh, W, c, meps, mmui, jsbow, e_exi, f1, bc);
[bc, W, e_exi, jsbow] = apply_bc(msh, bc, ebow2_bc, jsbow);
ebow2 = solve_helmholtz_2d_fd(msh, W, c, meps, mmui, jsbow, e_exi, f2, bc);

ebow = ebow1 + ebow2;
ebow_z = ebow((2/3)*length(ebow)+1:length(ebow));

%% Postprocessing
idx_edge_x = 1:msh.np;
idx_edge_y = 1+msh.np:2*msh.np;
idx_edge_z = 1+2*msh.np:3*msh.np;
ebow1_x = ebow1(idx_edge_x);
ebow1_y = ebow1(idx_edge_y);
ebow1_z = ebow1(idx_edge_z);
ebow1_abs = sqrt(abs(ebow1_x).^2 + abs(ebow1_y).^2 + abs(ebow1_z).^2);
ebow2_x = ebow2(idx_edge_x);
ebow2_y = ebow2(idx_edge_y);
ebow2_z = ebow2(idx_edge_z);
ebow2_abs = sqrt(abs(ebow2_x).^2 + abs(ebow2_y).^2 + abs(ebow2_z).^2);
ebow_x = ebow(idx_edge_x);
ebow_y = ebow(idx_edge_y);
ebow_z = ebow(idx_edge_z);
ebow_abs = sqrt(abs(ebow_x).^2 + abs(ebow_y).^2 + abs(ebow_z).^2);

idx = 1+NPML(2):length(ymesh)-NPML(4);  % Indices at which to evaluate the field
y = ymesh(idx);                             % y values at those indices

if plot_field
    figure
    [X,Y] = meshgrid(msh.xmesh, msh.ymesh);
    e_surf = reshape(ebow_abs, [msh.nx, msh.ny]);
    e_surf_plot = surf(X,Y,e_surf');
    %xlim([0, L])
    ylim([-h/2, h/2])
    set(e_surf_plot,'LineStyle','none')
    set(gca,'ColorScale','log')
    title('Absolute value of magnetic field','Interpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex')
    ylabel('$y$ (m)','Interpreter','latex')
    zlabel('Absolute value','Interpreter','latex')
end

% Intensity calculation % TODO: CAN'T add intensities!!!
idx_screen = msh.nx * (1:msh.ny) - NPML(1);
e1_screen = ebow1_abs(idx_screen)';
e2_screen = ebow2_abs(idx_screen)';
%e_screen = ebow(msh.nx * (1:msh.ny) - NPML(1))';
e1_screen = e1_screen(idx);
e2_screen = e2_screen(idx);
%e_screen = e_screen(1+NPML(2):end-NPML(4));
I1 = c0*eps0/2 * abs(e1_screen).^2;
I2 = c0*eps0/2 * abs(e2_screen).^2;
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

I1_farfield = intensity_farfield(E1, lambda1, delta, L, y);
I2_farfield = intensity_farfield(E2, lambda2, delta, L, y);
I_farfield = I1_farfield + I2_farfield;
I1_helmholtz = intensity_helmholtz(E1, lambda1, delta, L, y, ceil(length(idx_bc)));
I2_helmholtz = intensity_helmholtz(E2, lambda2, delta, L, y, ceil(length(idx_bc)));
I_helmholtz = I1_helmholtz + I2_helmholtz;
if plot_intensity_ana % TODO: FIX normalization
    plot(y, I_helmholtz/max(I_helmholtz), 'b--', 'DisplayName', 'Analytical (Helmholtz)')
    plot(y, I_farfield/max(I_farfield), 'r--', 'DisplayName', 'Analytical (farfield)')
end

% Error calculation
I_err = norm(I/max(I) - I_farfield/max(I_farfield))/norm(I_farfield/max(I_farfield));
I_err_helmholtz = norm(I/max(I) - I_helmholtz/max(I_helmholtz))/norm(I_helmholtz/max(I_helmholtz));
if calc_intensity_err
    fprintf('Relative L2 error between numerical and farfield solution = %f%% \n', 100*I_err)
    fprintf('Relative L2 error between numerical and Helmholtz solution = %f%%', 100*I_err_helmholtz)
end


