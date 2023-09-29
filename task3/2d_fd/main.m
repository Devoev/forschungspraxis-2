%% Task3: Thin-film interference

clc
clear all

% Path mesh functions
path_msh_func = '../fit/2d/mesh';
path_mat_func = '../fit/2d/matrices';
path_solver_func = '../fit/2d/solver';
path_util_func = '../fit/2d/util';
%path_verify_func = './task1/verifications';

% Add paths
cd('../');
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thin_film = 0;
plot_field = 1;
plot_intensity=1;     % c)


%% Problem Definition
c = 3e8;            % [m/s]
eps0 = 8.854e-12;
mui = 1/(4*pi*1e-7);
n = 2;

lambda1 = 430e-9;   % [m]
f1 = c/lambda1;     % [Hz]
omega1 = 2*pi*f1;   % [1/s]
E1 = 250;           % [V/m]

lambda2 = 510e-9;   % [m]
f2 = c/lambda2;     % [Hz]
omega2 = 2*pi*f2;   % [1/s]
E2 = 500;           % [V/m]

%         thin film geometry           |   
%        -h/2    L2[PML]  h/2          |      b: middle index
%         |------ y ------|            |         
%        # | | | | | | | | #           |      d: distance between the
%  [PEC] # v v v v v v v v #   [PEC]   |         excitations in wave-
%  L3    #                 #   L1      |         length  
%  x     #-----------------#border1(L/2)|
%        #----thin film----#  a        |      L: side index bc
% border2#-----------------#           |         [L1, L2, L3, L4]
%        #                 #           |
%        ###################           |
%               L4 [PML]

%% geometry
h = 4e-6;   % [m]
L = 10e-6;  % [m]
a = 100e-9; % [m]

NPML = [20, 0, 20, 0];  % [L1, L2, L3, L4]; 0,1:=PMC
pml_space = 20; % spacing cells between pml layer and "real" model layers


%% mesh
elem_per_wavelength = 10;
dx = lambda1*(4*pml_space)/elem_per_wavelength;  % Extra space in +x direction
xmesh = linspace(0, L + dx, ceil((L + dx)/lambda1*elem_per_wavelength));
ymesh = linspace(-h/2, h/2, ceil( h/lambda1*elem_per_wavelength ));
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
    eps1 = (n^2)*eps;
    eps_vec = sparse(ones(msh.np, 1))*eps;
    % set all thin film elements to eps1
    for idx = border1_x_idx:border2_x_idx
        x_layer = all_elem(~(mod(all_elem, idx)));
        eps_vec(x_layer) = eps1;
    end
    eps = eps_vec
else
    % eps scalar -> homogenous isotrope meps
    eps = eps0;
end

%% excitation
indices = 1:msh.np;
idx = msh.np+indices(~(mod(indices,msh.nx)-2*pml_space));  % idx of all L2 boundary elements
% idx = round(msh.np*0.5)+round(msh.nx*0.5);
jsbow = sparse(3*msh.np, 1);
ebow1_bc = NaN(3*msh.np, 1);
ebow2_bc = NaN(3*msh.np, 1);
ebow1_bc(idx) = E1;
ebow2_bc(idx) = E2;
% Anregung erstmal in z-Richtung (wie Task1). entspricht Aufgabe d)

%% Create matrices
[C, ~, ~] = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);

meps = createMeps_2D(msh, ds, da, dat, ones(msh.np, 1), eps);
mmui = createMmui_2D(msh, ds, dst, da, ones(msh.np, 1), mui);

%% solve system
% Solve

ebow1 = frequency_domain_2D(msh, C, meps, mmui, 0, jsbow, idx, ebow1_bc(idx), omega1, NPML);
ebow2 = frequency_domain_2D(msh, C, meps, mmui, 0, jsbow, idx, ebow2_bc(idx), omega2, NPML);
ebow = ebow1 + ebow2;


%% Postprocessing
if plot_field
    figure
    [X,Y] = meshgrid(msh.xmesh, msh.ymesh);
    e_surf = reshape(real(ebow(msh.np+1:2*msh.np)), [msh.nx, msh.ny]);
    e_surf_plot = surf(X,Y,e_surf');
    xlabel(' X ');
    ylabel(' Y ');
    %zlim([-E2*6 E2*6])
    set(e_surf_plot,'LineStyle','none')
    set(gca,'ColorScale','log')
end

% Intensity calculations % TODO: CAN'T add intensities!!!
xL_idx = all_elem(~(mod(all_elem, msh.nx)-msh.nx+2*pml_space));
x0_idx = all_elem(~(mod(all_elem, msh.nx)-2*pml_space));
%e2_screen = ebow2(msh.nx * (1:msh.ny) - offset(1))';
%e_screen = ebow(msh.nx * (1:msh.ny) - NPML(1))';
e_x0 = ebow(x0_idx);
e_xL = ebow(xL_idx);
% iam using scalar eps0 because both intensities are calculated outside of 
% the thin film medium
I_0 = c*eps0/2 * abs(e_x0).^2;
I_L = c*eps0/2 * abs(e_xL).^2;
%I = I1 + I2;

if plot_intensity
    figure
    plot(msh.ymesh, I_L/max(I_L), 'DisplayName', 'Numerical', 'color', '#1e8080')
    hold on
    title('Intensity of transmitted Wave at $x=L=10^6$m','Interpreter','latex')
    xlabel('Position at the screen $y$ (m)','Interpreter','latex')
    ylabel('Intensity $I$','Interpreter','latex')
    xlim([-h/2, h/2])
    legend()

    figure
    plot(msh.ymesh, I_0/max(I_0), 'DisplayName', 'Numerical', 'color', '#1e8080')
    hold on
    title('Intensity of reflected Wave at $x=0$m','Interpreter','latex')
    xlabel('Position at the screen $y$ (m)','Interpreter','latex')
    ylabel('Intensity $I$','Interpreter','latex')
    xlim([-h/2, h/2])
    legend()
end


