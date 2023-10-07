%% Task3: Thin-film interference

clc
clear all

filePath = matlab.desktop.editor.getActiveFilename;
baseDir = extractBefore(filePath,"task3");

% Paths to add
path_msh_func = append(baseDir, 'fit\2d\mesh');
path_mat_func = append(baseDir, 'fit\2d\matrices');
path_solver_func = append(baseDir, 'fit\2d\solver');
path_util_func = append(baseDir, 'fit\2d\util');
path_task = append(baseDir, 'task4\2d_fd');

% Add paths
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thin_film = 1;
plot_field = 1;
plot_analytic = 1;
plot_intensity = 0;     % c)
use_y_symmetry = 0; %ToDo add Symmetrie
polarisation = 'y'; 


%% Problem Definition
c = 3e8;            % [m/s]
eps0 = 8.854e-12;
mui0 = 1/(4*pi*1e-7);
n = 2;

lambda1 = 430e-9;   % [m]
f1 = c/lambda1;     % [Hz]
omega1 = 2*pi*f1;   % [1/s]
E1 = 250;           % [V/m]

lambda2 = 510e-9;   % [m]
f2 = c/lambda2;     % [Hz]
omega2 = 2*pi*f2;   % [1/s]
E2 = 500;           % [V/m]

%         thin film geometry                |   
%              L3[PEC]                      |      b: middle index
%                                           |         
%        #####################----y=-h/2    |      d: distance between the
%  [PML] ->      |   |       #   [PML]      |         excitations in wave-
%  L4    ->      |   |       #   L2         |         length  
%  x=0   ->      |   |       #   x=L        |
%        ->      |   |       #              |      L: side index bc
%        ->      |   |       #              |         [L1, L2, L3, L4]
%        ->      |   |       #              |
%        #####################---- y=h/2    |
%               L/2  L/2s+a    
%              L1 [PEC]

%% geometry
h = 4e-6;   % [m]
L = 10e-6;  % [m]
a = 100e-9; % [m]

NPML = [20, 20, 20, 20];  % [L1, L2, L3, L4]; 0,1:=PMC
pml_space = 20; % spacing cells between pml layer and "real" model layers
bc.bc = ["PEC", "OPEN", "PEC", "OPEN"];   % Total offset from boundaries
bc.NPML = NPML;


%% mesh
elem_per_wavelength = 15;
dx = lambda1*(NPML(2) + NPML(4))/elem_per_wavelength;  % Extra space in +x direction
xmesh = linspace(0, L + dx, ceil((L + dx)/lambda1*elem_per_wavelength));
ymesh = linspace(-h/2, h/2, ceil( h/lambda1*elem_per_wavelength )); % No Extra space in y dir. (PEC)
msh = cartMesh_2D(xmesh, ymesh);

% border x indices of different permittivity == thin film borders
border1_x_idx = round(msh.nx*0.5); 
border2_x_idx = border1_x_idx + round(a*msh.nx/(L+dx));
actual_thickness = (border2_x_idx-border1_x_idx)*((L+dx)/msh.nx)

%% thin film material
% refraction idx n=2: wavevelocity in thin film 0.5*c0 
%  -> multiplying mu or eps value for thin film area by 4
all_elem = 1:msh.np;
if thin_film
    % create inhomogenous eps vector
    eps1 = (n^2);
    eps_vec = sparse(ones(msh.np, 1))*1;
    % set all thin film elements to eps1
    for idx = border1_x_idx:border2_x_idx
        %x_layer = all_elem(~(mod(all_elem, idx)));
        x_layer = idx + (0:(msh.ny-1))*(msh.nx);
        eps_vec(x_layer) = eps1;
    end
    eps = eps_vec;
else
    % eps scalar -> homogenous isotrope meps
    eps = sparse(ones(msh.np, 1))*1;
end

%% excitation

%Excitation in wohle distance h
idx_bc = calc_slit_idx(msh, h, use_y_symmetry, polarisation) + NPML(3); % Transform y-indices to canonical index
jsbow = sparse(3*msh.np, 1);
ebow1_bc = NaN(3*msh.np, 1);
ebow2_bc = NaN(3*msh.np, 1);
ebow1_bc(idx_bc) = E1;
ebow2_bc(idx_bc) = E2;

%% Create matrices

[C, ~, ~] = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);

%% ToDo use Box mesher 

meps = createMeps_2D(msh, ds, da, dat, eps, eps0);
mmui = createMmui_2D(msh, ds, dst, da, ones(msh.np, 1), mui0);

%% solve system
% Solve

[bc, W, e_exi, jsbow] = apply_bc(msh, bc, ebow1_bc, jsbow);
ebow1 = solve_helmholtz_2d_fd(msh, W, C, meps, mmui, jsbow, e_exi, f1, bc);
[bc, W, e_exi, jsbow] = apply_bc(msh, bc, ebow2_bc, jsbow);
ebow2 = solve_helmholtz_2d_fd(msh, W, C, meps, mmui, jsbow, e_exi, f2, bc);
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
ebow_abs = abs(abs(ebow_x).^2 + abs(ebow_y).^2 + abs(ebow_z).^2);

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
    title('Absolute value of electric field','Interpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex')
    ylabel('$y$ (m)','Interpreter','latex')
    zlabel('Absolute value','Interpreter','latex')
end

% Intensity calculations % TODO: CAN'T add intensities!!!

%Indices for the screen
xL_idx = msh.nx * (1:msh.ny) - NPML(1);
x0_idx = (1:msh.ny) + NPML(3);
e_x0 = ebow_abs(x0_idx);
e_xL = ebow_abs(xL_idx);
e_x0 = e_x0(idx);
e_xL = e_xL(idx);
% iam using scalar eps0 because both intensities are calculated outside of 
% the thin film medium
I_0 = c*eps0/2 * abs(e_x0).^2;
I_L = c*eps0/2 * abs(e_xL).^2;
%I = I1 + I2;

if plot_intensity
    figure
    plot(y, I_L, 'DisplayName', 'Numerical', 'color', '#1e8080')
    hold on
    title('Intensity of transmitted Wave at $x=L=10^6$m','Interpreter','latex')
    xlabel('Position at the screen $y$ (m)','Interpreter','latex')
    ylabel('Intensity $I$','Interpreter','latex')
    xlim([-h/2, h/2])
    legend()

    figure
    plot(y, I_0, 'DisplayName', 'Numerical', 'color', '#1e8080')
    hold on
    title('Intensity of reflected Wave at $x=0$m','Interpreter','latex')
    xlabel('Position at the screen $y$ (m)','Interpreter','latex')
    ylabel('Intensity $I$','Interpreter','latex')
    xlim([-h/2, h/2])
    legend()
end

if plot_analytic
    e_analytic = analytic_sol(msh.xmesh, E1, E2, lambda1, lambda2, L/2, actual_thickness, n);

    figure
    plot(msh.xmesh, abs(e_analytic),'DisplayName', 'Numerical', 'color', '#1e8080')
    hold on
    title('Absolute E-Field across x-Axis (analytic)','Interpreter','latex')
    xlabel('x Position (m)','Interpreter','latex')
    ylabel('absolut e-Field $|E|$','Interpreter','latex')
    xlim([0, 1.2*L])
    

end


