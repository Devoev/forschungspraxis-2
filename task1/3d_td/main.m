% Main 3d Solver
clc 
clear all

% Path mesh functions
path_fit_func = '../fit_3d';

% Add paths
cd('../');
addpath(path_fit_func)

%% Here are some options
bcalcMaxTimestep    = 0;
b2dFieldPlot        = 0;
bPlotonScreen       = 1;
bPlotScreen_mean    = 1;

%% Problem definition

c0 = 3e8;            % m/s
eps = 8.854e-12;
mu = 4e-7*pi;
mui = 1/(4*pi*1e-7);

lambda1 = 430e-9;   % m
f1 = c0/lambda1;     % Hz
omega1 = 2*pi*f1;   % 1/s
E1 = 250;           % V/m

%lambda2 = 510e-9;
%f2 = c/lambda2;
%omega2 = 2*pi*f2;
%E2 = 500;

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

d = 4e-6;       % slit distance
delta = 1e-6;   % slit width
h = 8e-6;       % domain in y dir. (parallel to source and screen)
L = 10e-6;      % domain in x dir. (direction source to screen)

% set open boundary for each side (overrides the pmc/pec bc)
open_bc = [true, true, true, true];  % [L1, L2, L3, L4];


%% Generate Mesh
elem_per_wavelength = 15;

%x_mesh = [linspace(-wavelength_offset_x * elem_per_wavelength * edge_size, 0, wavelength_offset_x * elem_per_wavelength), ...
%    linspace(0, L, ceil(L/edge_size))];

% x-dir: subtracted from distances (source, screen) (ToDo change that it is
% added)(ToDo change distance and time integration)
wavelength_offset_x = 2;      % defines distance from screen and source to boundary in x dir
elem_offset_x = wavelength_offset_x * elem_per_wavelength;  % distance in elements in x

% y-dir: added to domain wavelengthoffset
wavelength_offset_y = 10;     % defines offset from screen and source to boundary in (+/-) y dir
off_y = wavelength_offset_y * lambda1;
elem_offset_y = wavelength_offset_y * elem_per_wavelength; % distance in elements in y

  nz = 2;
  xmesh = linspace(0, L, (L/lambda1)*elem_per_wavelength);
  ymesh = linspace(-h/2 - off_y, h/2 + off_y, ((h+2*off_y)/lambda1)*elem_per_wavelength);
  zmesh = linspace(0,1e-7,nz);  % ToDo - set value for distance in z dir.
  msh = cartMesh(xmesh, ymesh, zmesh); 

  nx = msh.nx;
  ny = msh.ny;
  nz = msh.nz;

% grid distances (equidistant)
  delta_x = L/(nx-1);
  delta_y = (h+2*off_y)/(ny-1);
  delta_z = 1e-7/(nz-1);

Mx = msh.Mx;
My = msh.My;
Mz = msh.Mz;
np = msh.np;

%% Generate topological und geometric matrices
[c, s, st] = createTopMats(msh);
[ds, dst, da, dat] = createGeoMats(msh);

bcs = [ 1, 1, 1, 1, 1, 1];

Mmui = createMmui(msh, ds, dst, da, mui, bcs);
Mmu = nullInv(Mmui);

Meps = createMeps(msh, ds, da, dat, eps, bcs);
Mepsi = nullInv(Meps);

% init mur edges
[mur_edges, mur_n_edges, mur_deltas] = initMur(msh, open_bc);
%% CFL-condition for max timestep

 deltaTmaxCFL = sqrt(mu*eps)*sqrt(1/((1/delta_x^2)+(1/delta_y^2)+(1/delta_z^2)));
 fprintf('According to CFL: deltaTmax = %e\n',deltaTmaxCFL);

%% Max timestep via eigenvalue decomposition

if bcalcMaxTimestep
  
    A12 = -Mmui*c;
    A21 = (nullInv(Meps))*transpose(c);
    zero = sparse(3*np, 3*np);
    A = [zero, A12; A21, zero];
    [~, lambdaMaxA] = eigs(A,1);
    lambdaMaxA = abs(lambdaMaxA);
    % Workaround for octave
    if ~isscalar(lambdaMaxA)
        lambdaMaxA = lambdaMaxA(1,1);
    end

 deltaTmaxEigA = 2/lambdaMaxA;
 fprintf('From Eigenvalue calculation of A: deltaTmax = %e\n', deltaTmaxEigA);

end

%% Simulation parameter

 dt = 0.5*deltaTmaxCFL;
% sigma = 25*deltaTmaxCFL
% tend = 100*sigma;
 tend = 3.33e-14;  % t = L/c
 steps = ceil(tend/dt);
 sourcetype= 2;  % 1: Gauss Anregung, 2: Harmonisch, 3: Konstante Anregung

%% Excitation with current desnity

% Calculate BC indices
y_slit = [(-d - delta)/2, (-d + delta)/2, (d - delta)/2, (d + delta)/2]; % y values of upper and lower slit.
for i = 1:length(y_slit)
    % Find y-index closest to actual y_slit value
    [~,y_idx(i)] = min(abs(msh.ymesh - y_slit(i)));
end
y_idx = [y_idx(1):y_idx(2), y_idx(3):y_idx(4)]; % Find all y-indices between slits

% Set rhs and bc vectors
idx = elem_offset_x + msh.nx * (y_idx-1) + 2*np; % Transform y-indices to canonical index

%indices for excitation are just provisional
jsbow_space = zeros(3*np, 1);
jsbow_space(idx) = 1;
    %jsbow_space(2*np+ceil(0.5*nx)*Mx+ceil(0.5*ny)*My) = 1;
%    jsbow_space(2*np+2*Mx+2*My) = 1;
%    jsbow_space(2*np+nx-2*Mx+2*My) = 1;

% Gauss 
 jsbow_gauss = @(t)(jsbow_space * exp(-4*(((t-sigma)/sigma)^2)));

% Harmonic
% jsbow_harm = @(t)(jsbow_space * sin(pi*(t/sigma)));
 jsbow_harm = @(t)(jsbow_space * sin(omega1*t));

% Const
 jsbow_const = @(t)(jsbow_space * 1);

%% Initializations
ebow_new = sparse(3*np,1);
hbow_new = sparse(3*np,1);

% initialize figures
if b2dFieldPlot
    figure(1)
end

if bPlotonScreen
    figure(2)
end

draw_only_every = 4;


%% Calculation (and parts of postprocessing) for every timestep
t=0;
integration_step = 0;

% Time integration
for ii = 1:steps

    % calc time
      t = t + dt;

    % save old values
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    if sourcetype == 1
        % get excitation current
        if t <= 2*sigma
            js = jsbow_const(t);
        else
            js = zeros(3*np, 1);
        end
    elseif sourcetype == 2
        % harmonic excitation
        js = jsbow_harm(t);
    elseif sourcetype == 3
        % const excitation
        js = jsbow_const(t);
    end
    
    % do the leapfrog step
    [hbow_new,ebow_new] = leapfrog(hbow_old, ebow_old, js, Mmui, Mepsi, c, dt);
    % apply open boundary w. mur cond
    ebow_new = applyMur(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt);

    % calculate intensity out of e_bow
    I = c0*eps/2 * abs(ebow_new).^2;
    
    % plots
    if mod(ii, draw_only_every)

        if b2dFieldPlot             %show e-field on x,y surface
		    e_surface = reshape(ebow_new((2*np+1):(2*np+1*Mz)),nx,ny);
            figure(1)
                mesh(e_surface)
                axis([1 nx 1 ny -400 400])
                caxis([-100 100])
                title('e-bow on x-y surface','Interpreter','latex')
                xlabel('x','Interpreter','latex')
                ylabel('y','Interpreter','latex')
                zlabel('$e-bow$','Interpreter','latex')
        end 
       
        if bPlotonScreen        % show e-field on the screen
                %e_screen = ebow_new((2*np+nx*(ny-5)+1):(2*np+nx*(ny-4)))';
                % screen is elem_offset_x ahead of bc
                %e_screen = ebow_new(-elem_offset_x + msh.nx * (1:msh.ny) + 2*np);
                I_screen = I(-elem_offset_x + msh.nx * (elem_offset_y:(msh.ny-elem_offset_y)) + 2*np);
                figure(2)
                    plot(ymesh(elem_offset_y:(length(ymesh)-elem_offset_y)),I_screen)
                    ylim([0 c0*eps/2*E1^2])
                    xlim([ymesh(1)+off_y ymesh(end)-off_y])
                    title('Intensity at the screen at $x=L=10^6$m','Interpreter','latex')
                    xlabel('Position at the screen $y$ (m)','Interpreter','latex')
                    ylabel('$I$','Interpreter','latex')
        end

        drawnow
     
        % TODO: Calculate time avareged intensity (before reflections come in)
        % time when wave arrives at screen, source and screen are 1 wavelengthoffset away from bc
        % wavelengthoffset need to be bigger than 1 such that there are no
        % reflections from the back
        time_arrive = (L - 2*wavelength_offset_x*lambda1)/c0;
        % time when 1st full wave passed screen
        time_pass = time_arrive + 1/f1;
        if (time_arrive < t) & (t < time_pass)
            integration_step = integration_step + 1;
            I_screen = I(-elem_offset_x + msh.nx * (elem_offset_y:(msh.ny-elem_offset_y)) + 2*np);
            I_screen_sample(integration_step,:) = I_screen;
        end
    end
end

if bPlotScreen_mean
    I_screen_mean = mean(I_screen_sample);
    figure(3)
    plot(ymesh(elem_offset_y:(length(ymesh)-elem_offset_y)),I_screen_mean)
    %ylim([0 c0*eps/2*E1^2])
    xlim([ymesh(1)+off_y ymesh(end)-off_y])
    title('I mean at the screen at $x=L=10^6$m','Interpreter','latex')
    xlabel('Position at the screen $y$ (m)','Interpreter','latex')
    ylabel('$I$','Interpreter','latex')
end
