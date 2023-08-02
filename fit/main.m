% Path mesh functions
path_msh_func = '.\mesh';
path_mat_func = '.\matrices';
path_solver_func = '.\solver';
path_util_func = '.\util';
path_verify_func = '.\verifications';

% Add paths
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func, path_verify_func)

%% Problem Definition
f = 50e6;
omega = 2*pi*f;
excitation = 1000;
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

elem_per_wavelength = 20;

x = 80;
y = 100;
d = 20; % "Gitter-Abstand"

lambda = 3e8/f ; % Wavelength;
[ msh ] = cartMesh2D( linspace(1,x*lambda, x*elem_per_wavelength), linspace(1,y*lambda, y*elem_per_wavelength) );

b = floor((x*elem_per_wavelength)/2) + msh.ny * (20 + elem_per_wavelength);
jsbow = sparse(msh.np, 1);
jsbow(floor(b - d*elem_per_wavelength ),1) = excitation;
jsbow(ceil(b + d*elem_per_wavelength),1) = excitation;

display(fresnel_number(d, y*elem_per_wavelength, max(diff(msh.xmesh)), max(diff(msh.ymesh)), f))

%% Solution
eps = 8.854e-12;
mui = 1/(4*pi*1e-7);

ebow = solveHelmholtzTE(msh, eps, mui, jsbow, omega, 0);
save('ebow.mat', 'ebow')
ibov = reshape(real(ebow*exp(-1i*omega)),[msh.nx, msh.ny]);

%% Postprocessing
figure
[X,Y] = meshgrid(msh.xmesh, msh.ymesh);
h = surf(X,Y,ibov');
set(h,'LineStyle','none')
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

E_ana = helmholtz_analytic(lambda, x, y, elem_per_wavelength, omega, excitation);
figure 
plot(1:length(E_ana), abs(E_ana))

%%
function fnum = fresnel_number(nD ,nL, dx, dy, f)
    a = nD * dx;
    L = nL * dy;
    fnum = a^2 * f / (L * 3e8);
    assert(fnum < .2);
end
