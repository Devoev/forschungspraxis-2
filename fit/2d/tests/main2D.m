%% Clear workspace

clc
clear

f = 2e9;

%% Define calculation domain

% Define mesh coordinates in x- and y-direction
xmesh = linspace(0,1,61);
ymesh = linspace(0,1,61);


% Define relative permittivity and inverse relative permeability
epsilon = 1;
mui = 1;


%% Create Mesh

msh = cartMesh(xmesh, ymesh);

nx = msh.nx;
ny = msh.ny;
np = msh.np;

[ds, dst, da, dat] = createGeoMats(msh);

idx_dof = getNotGhostEdges(msh);

meps = createMeps(msh, ds, dat, epsilon);

mkap = createMeps(msh, ds, dat, 4000000000);


mmui = createMmui(msh, dst, da, mui);

[c, s, st] = createTopMats(msh);

jsbow = sparse(3*msh.np, 1);
jsbow(9365) = 1;
omega = 2 * pi * f;
A = -c'*mmui*c + omega^2*meps - 1i*omega*mkap;
b = 1j*omega*jsbow;

% Deflate system matrix

idx_dof = nonzeros(idx_dof);

b = b(idx_dof);
A = A(idx_dof, idx_dof);

% solve equation
%    [ebow, flag, relRes, iter, resVec] = gmres(A, rhs, 20, 1e-10, 1000); % TODO: direct vs iteratve?
%    if flag == 0
%      fprintf('gmres(20): converged at iteration %2d_te to a solution with relative residual %d.\n',iter,relRes);
%    else
%      error('gmres(20): some error ocurred, please check flag output.')
%    end
%    relRes = resVec./norm(rhs);
ebow = sparse(3*msh.np, 1);
ebow(idx_dof) = A\b;

ebow = real(ebow);

zlimit = 700;

z_plane = 1;
idx2plot = 2*np+1:3*np;
ebow_mat = reshape(ebow(idx2plot),nx,ny);
figure(1)
mesh(ebow_mat)
xlabel('i')
ylabel('j')
zlabel(['z-Komponente des E-Feldes für z=',num2str(z_plane)])
axis([1 nx 1 ny -zlimit zlimit])

drawnow



