% Path mesh functions
path_msh_func = '..\mesh';

path_mat_func = '..\matrices'

f = 1e6;



% Add paths
addpath(path_msh_func, path_mat_func)


 elm = 100;

[ msh ] = cartMesh2D( linspace(1,elm,20*elm), linspace(1,elm,20*elm) );

jsbow = sparse(msh.np, 1);
% jsbow(2000*2000/2 - 1000,1) = 1000;
jsbow(1,1) = 1000;
jsbow(20*elm,1) = 1000;


% Anzahl der Rechenpunkte des Gitters
np = msh.np;

[c, g, st] = createTopMats2DTE(msh);


% TODO: 2D material matrices
meps = 8.854e-12 * speye(msh.np, msh.np);
mmui = 1/(4*pi*1e-7) * speye(2*msh.np, 2*msh.np);

% Berechnung der Kreisfrequenz
omega = 2*pi*f;

% Berechnung Systemmatrix A und rechte Seite rhs
A = -c'*mmui*c + omega^2*meps;
A = st * mmui * g + omega^2*meps;
rhs = 1j*omega*jsbow;

% solve equation
ebow = A\rhs; % TODO: direct vs iteratve?


% Post processing
bbow = -c*ebow / (1i*omega);
hbow = mmui*bbow;

ibov = reshape(real(ebow*exp(-1i*omega)),[msh.nx, msh.ny]);


%spy(ibov)

imagesc(ibov)
colorbar



