% Path mesh functions
path_msh_func = '..\mesh';

path_mat_func = '..\matrices';
path_solver_func = '..\solver';

% Add paths
addpath(path_msh_func, path_mat_func, path_solver_func)

f = 50e6;
omega = 2*pi*f;
elm = 100;

[ msh ] = cartMesh2D( linspace(1,elm,20*elm), linspace(1,elm,20*elm) );

eps = 8.854e-12;
mui = 1/(4*pi*1e-7);

d = 200; % "Gitter-Abstand"
jsbow = sparse(msh.np, 1);
% jsbow(2000*2000/2 - 1000,1) = 1000;
jsbow(10*elm - d,1) = 1000;
jsbow(10*elm + d,1) = 1000;

ebow = solveHelmholtz2D(msh, eps, mui, jsbow, omega, 0);

ibov = reshape(real(ebow*exp(-1i*omega)),[msh.nx, msh.ny]);


%spy(ibov)

imagesc(ibov)
colorbar

figure

[X,Y] = meshgrid(msh.xmesh, msh.ymesh);
h = surf(X,Y,ibov);
set(h,'LineStyle','none')

figure

intensity = ibov(:,end);
plot(1:length(intensity), intensity)
