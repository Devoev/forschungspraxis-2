% Path mesh functions
path_msh_func = './fit/mesh';
path_mat_func = './fit/matrices';
path_solver_func = './fit/solver';
path_util_func = './fit/util';
path_verify_func = './task1/verifications';

% Add paths
cd('../');
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func, path_verify_func)

%% Problem Definition
c = 3e8;            % m/s

lambda1 = 430e-9;   % m
f1 = c/lambda1;     % Hz
omega1 = 2*pi*f1;   % 1/s
E1 = 250;           % V/m

lambda2 = 510e-9;
f2 = c/lambda2;
omega2 = 2*pi*f2;
E2 = 500;

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

delta = 1e-6;   % slit width
d = 4e-6;       % slit distance
h = 8e-6;       % screen height
L = 10e-6;      % screen distance


%% Generate Mesh
elem_per_wavelength = 10;
xmesh = linspace(0, L, L/lambda1*elem_per_wavelength);
ymesh = linspace(-h/2, h/2, h/lambda1*elem_per_wavelength);
msh = cartMesh2D(xmesh, ymesh);

[~,y1] = min(abs(msh.ymesh - d/2)); % Index of 1st slit % TODO: Increase slit width by 'delta'
[~,y2] = min(abs(msh.ymesh + d/2)); % Index of 2nd slit

jsbow = sparse(msh.np, 1);
idx1 = msh.nx * (y1-1); %  TODO: Add offset for PML boundary
idx2 = msh.nx * (y2-1);
jsbow(idx1) = E1; % TODO: Use J instead of E
jsbow(idx2) = E1;


fprintf('Fresnel number = %f', fresnel_number(delta, L, lambda1))

[X,Y] = meshgrid(xmesh, ymesh);
plot(X, Y, 'blue', X', Y', 'blue')
hold on
plot(0, d/2, 'x', 'color', '#A2142F', 'linewidth', 6)
plot(0, -d/2, 'x', 'color', '#A2142F', 'linewidth', 6)
xline(0, 'red')
xline(L, 'red')
yline(-h/2, 'red')
yline(h/2, 'red')


%% Solution
eps = 8.854e-12;
mui = 1/(4*pi*1e-7);

ebow = solveHelmholtzTE(msh, eps, mui, jsbow, omega1, 0);
save('ebow.mat', 'ebow')
ibov = reshape(real(ebow*exp(-1i*omega1)),[msh.nx, msh.ny]);

%% Postprocessing
figure
[X,Y] = meshgrid(msh.xmesh, msh.ymesh);
hbow = surf(X,Y,ibov');
set(hbow,'LineStyle','none')
set(gca,'ColorScale','log')
figure

% TODO: proper intensity calculation
intensity = ibov(20:end-20,end-20).^2;
plot(1:length(intensity), abs(intensity))



%% verifications

%analytical solution of Helmholtz eq:
%for easy comparison: following excitation is assumed:
%jsbow(floor(b),1) = excitation; (1 excitation in the center)
%following "screen" is used
%intensity = ibov(20:end-20,end-20).^2;
%formula is described in LaTEx

E_ana = helmholtz_analytic(lambda1, L, h, elem_per_wavelength, omega1, E1);
figure 
plot(1:length(E_ana), abs(E_ana))
