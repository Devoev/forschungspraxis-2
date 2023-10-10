%% Add paths

% Clear variables
clc
clear

% Get parent directory
filePath = matlab.desktop.editor.getActiveFilename;
[ParentFolderPath] = fileparts(filePath);
parent = fileparts(ParentFolderPath) ;

% Paths to add
path_msh_func = append(parent, '\mesh');
path_mat_func = append(parent, '\matrices');
path_solver_func = append(parent, '\solver');
path_util_func = append(parent, '\util');

% Add paths
addpath(path_msh_func, path_mat_func, path_solver_func, path_util_func)


%% Edit basic calculation domain 

% Steps of mesh
xmesh = linspace(-0.5,3.5,404);
ymesh = linspace(-0.5,0.5,101);

% Create basic mesh object
msh = cartMesh_2D(xmesh, ymesh); 
Mx = msh.Mx;
My = msh.My;
nx = msh.nx;
ny = msh.ny;
np = msh.np;

% Edit boundary conditions for polarization of E in z-direction
bc.bc = ["PMC", "OPEN", "PMC", "OPEN"];


%% Edit material regions and add them to the object material_regions

% Add basic constants to material_regions object
material_regions.epsilon0 = 8.854187e-12;
material_regions.mu0i = 1/(pi*4e-7);

% Regions for relative permittivity
% Relative permittivity everywhere equal to one
boxesEpsilonR(1).box = [1, nx, 1, ny];
boxesEpsilonR(1).value = 1;
material_regions.boxesEpsilonR = boxesEpsilonR;

% Regions for inverse relative permeability
% Inverse relative permeability everywhere equal to one
boxesMuiR(1).box = [1, nx, 1, ny];
boxesMuiR(1).value = 1;
material_regions.boxesMuiR = boxesMuiR;


%% Create Excitation for TEM wave with two polarizations in y and z

% Edge size in y-direction and x-direction
dx = xmesh(2) - xmesh(1);
dy = ymesh(2) - ymesh(1);

% Frequency for harmonic excitation
f = 1e9;

% Function for harmonic excitation in time domain
e_harm = @(t)(sin(2*pi*f*t));

% Time step size for harmonic excitation
dt = 1e-11;

% End time for harmonic excitation
t_end = 6.25/f;

% Create empty jsbow_space vector
jsbow_excitation = NaN(3*np, 1);

% Create excitation vector for the electric field
e_exitation = NaN(3*np, 1);

% Determine indices n for electric field excitation
x_h = ceil(max(size(xmesh))/2)+1;
n = 1 + (x_h-1)*Mx + ((1:ny)-1)*My; 

% Set corresponding entries in e_exi to the value of the desired amplitude
% -> Use corresponding edges as source
% Polarization of E in z-direction
e_exitation(n+2*np) = 1;  
% Polarization of E in y-direction
% e_exitation(n+1*np) = 1*dy;  


%% Apply boundary conditions and get excitation vectors for the simulation

[bc, W, e_exi, jsbow] = apply_bc(msh, bc, e_exitation, jsbow_excitation);


%% Generate matrices for calculation

[MAT, bc] = generate_MAT(msh, bc, material_regions, ["CurlP"]); %#ok<NBRAK2> 


%% Simulate in frequency domain

[ebow_freq, hbow_freq] = solve_helmholtz_2d_fd(msh, W, MAT.c, MAT.meps, MAT.mmui, jsbow, e_exi, f, bc);


%% Calculate power emitted through each surface

S_freq = CalcPowerSurfaceXY(msh, ebow_freq, hbow_freq, MAT.ds, MAT.dst, MAT.da);


%% Simulate in time domain

% Initialize open boundary condition if needed
[mur_edges,mur_n_edges, mur_deltas] = initMur_2D(msh, bc);

% Initialize ebow and hbow
ebow_new = sparse(3*np,1);
hbow_new = sparse(3*np,1);

% Add inverse permittivity matrix
MAT.mepsi = nullInv(MAT.meps);

% Initialize calculation of avarage power
S_td = zeros(3*np,1);
i_steps = 0;

% Calculate time steps
for t = linspace(0,t_end,ceil(t_end/dt))

    % Old values
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    % Calculate value for excitation
    e_exi_old = e_exi * e_harm(t);
    e_exi_new = e_exi * e_harm(t+1);

    % Execute timestep with leapfrog
    [ebow_new,hbow_new] = solve_leapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.c,dt,W);

    % Apply open boundary with mur condition
    ebow_new = applyMur_2D(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt, bc);

    % Sum up power through each surface
    if t >= t_end - 1/f
        S_td = S_td + CalcPowerSurfaceXY(msh, ebow_new, hbow_new, MAT.ds, MAT.dst, MAT.da);
        i_steps = i_steps + 1;
    end

end

% Get time avarage of the power through each surface
S_td = S_td/i_steps;


%% Plot results for wave with polarization of E in z-direction 

% Indices for edges of quantites to plot
x_indices = (51:max(size(xmesh))-50);
y_indices = ceil(max(size(ymesh))/3);
idx_plot =  1 + (x_indices-1) * Mx + (y_indices-1) * My;

% Analytical solution for electric and magnetic field field
k = 2*pi*f*sqrt(MAT.epsilon0/MAT.mu0i);
Z = sqrt(1/MAT.epsilon0/MAT.mu0i);
x_offste = xmesh(x_h);
e_analytic = cos(k * (xmesh(x_indices)-x_offste));
h_analytic = cos(k * (xmesh(x_indices)-x_offste+dx/2))/Z;
h_analytic(x_h-51+1:end) = -h_analytic(x_h-51+1:end);

% Analytical solution for power emitted in negative x-direction
power_neg_x = -0.5 / Z * (ny-3)*dy;

% Analytical solution for power emitted in positive x-direction
power_pos_x = 0.5 / Z * (ny-3)*dy;

% Plot electric field
figure(1)
subplot(2,1,1);
plot(xmesh(x_indices), e_analytic, xmesh(x_indices), real(ebow_freq(idx_plot+2*np)), xmesh(x_indices), ebow_new(idx_plot+2*np));
xlim([0  3]);
ylim([-1.5 1.5]);
legend({'Analytic','Numerical FD', 'Numerical TD'},'Location','southwest');
title('E polarized in z-direction: E in z-direction')

% Plot magnetic field
subplot(2,1,2);
plot(xmesh(x_indices)-dx/2,h_analytic,xmesh(x_indices)-dx/2,real(hbow_freq(idx_plot+np))/dy,xmesh(x_indices)-dx/2,hbow_old(idx_plot+np)/dy);
xlim([0  3]);
ylim([-1.5/Z  1.5/Z]);
legend({'Analytic','Numerical FD', 'Numerical TD'},'Location','southwest');
title('E polarized in z-direction: H in y-direction')
drawnow


%% Calculate errors for wave with polarization of E in z-direction 

% Errors electric field
error_E_FD = norm(real(ebow_freq(idx_plot+2*np)) - e_analytic') / norm(e_analytic);
error_E_TD = norm(ebow_new(idx_plot+2*np) - e_analytic') / norm(e_analytic);

% Errors magnetic field
error_H_FD = norm(real(hbow_freq(idx_plot+np))/dy - h_analytic') / norm(h_analytic);
error_H_TD = norm(hbow_old(idx_plot+np)/dy - h_analytic') / norm(h_analytic);

% Display errors:
disp('Errors for electromagnetic field with polarization of E in z-direction:');
disp(['Error electric field FD: ', num2str(error_E_FD)]);
disp(['Error electric field TD: ', num2str(error_E_TD)]);
disp(['Error magnetic field FD: ', num2str(error_H_FD)]);
disp(['Error magnetic field TD: ', num2str(error_H_TD)]);


%% Compare numerical calculated emitted power with analytical results

% Calculate error power emitted to yz-plane at xmesh(52) 
n_power = 1 + (52-1) * Mx + ((1:ny)-1) * My + np;
relative_err_neg_x_fd = abs((sum(real(S_freq(n_power))) - power_neg_x))/abs(power_neg_x);
relative_err_neg_x_td = abs((sum(S_td(n_power)) - power_neg_x))/abs(power_neg_x);

% Calculate error power emitted to yz-plane at xmesh(354) 
n_power = 1 + (354-1) * Mx + ((1:ny)-1) * My + np;
relative_err_pos_x_fd = abs((sum(real(S_freq(n_power))) - power_pos_x))/abs(power_pos_x);
relative_err_pos_x_td = abs((sum(S_td(n_power)) - power_pos_x))/abs(power_pos_x);

% Display errors:
disp('Relative errors for emitted power in positive and negative x-direction:');
disp(['Relative error negative x-direction FD: ', num2str(relative_err_neg_x_fd)]);
disp(['Relative error negative x-direction TD: ', num2str(relative_err_neg_x_td)]);
disp(['Relative error positive x-direction FD: ', num2str(relative_err_pos_x_fd)]);
disp(['Relative error positive x-direction TD: ', num2str(relative_err_pos_x_td)]);


% %% Edit boundary conditions for polarization of E in y-direction
% bc.bc = ["PEC", "OPEN", "PEC", "OPEN"];
% 
% 
% %% Apply boundary conditions and get excitation vectors for the simulation
% 
% [bc, W, e_exi, jsbow] = apply_bc(msh, bc, e_exitation, jsbow_excitation);
% 
% 
% %% Generate matrices for calculation
% 
% [MAT, bc] = generate_MAT(msh, bc, material_regions, ["CurlP"]); %#ok<NBRAK2> 
% 
% 
% %% Simulate in frequency domain
% 
% [ebow_freq, hbow_freq] = solve_helmholtz_2d_fd(msh, W, MAT.c, MAT.meps, MAT.mmui, jsbow, e_exi, f, bc);
% 
% 
% %% Simulate in time domain
% 
% % Initialize open boundary condition if needed
% [mur_edges,mur_n_edges, mur_deltas] = initMur_2D(msh, bc);
% 
% % Initialize ebow and hbow
% ebow_new = sparse(3*np,1);
% hbow_new = sparse(3*np,1);
% 
% % Add inverse permittivity matrix
% MAT.mepsi = nullInv(MAT.meps);
% 
% % Plot parameter for "movie"
% 
% % Calculate time steps
% for t = linspace(0,t_end,ceil(t_end/dt))
% 
%     % Old values
%     ebow_old = ebow_new;
%     hbow_old = hbow_new;
% 
%     % Calculate value for excitation
%     e_exi_old = e_exi * e_harm(t);
%     e_exi_new = e_exi * e_harm(t+1);
% 
%     % Execute timestep with leapfrog
%     [ebow_new,hbow_new] = solve_leapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,MAT.mmui,MAT.mepsi,MAT.c,dt,W);
% 
%     % Apply open boundary with mur condition
%     ebow_new = applyMur_2D(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt, bc);
% 
% end
% 
% 
% %% Plot results for wave with polarization of E in y-direction 
% 
% % Indices for edges of quantites to plot
% x_indices = (51:max(size(xmesh))-50);
% y_indices = ceil(max(size(ymesh))/3);
% idx_plot =  1 + (x_indices-1) * Mx + (y_indices-1) * My;
% 
% % Analytical solution for electric and magnetic field field
% k = 2*pi*f*sqrt(MAT.epsilon0/MAT.mu0i);
% Z = sqrt(1/MAT.epsilon0/MAT.mu0i);
% x_offste = xmesh(x_h);
% e_analytic = cos(k * (xmesh(x_indices)-x_offste));
% h_analytic = -cos(k * (xmesh(x_indices)-x_offste+dx/2))/Z;
% h_analytic(x_h-51+1:end) = -h_analytic(x_h-51+1:end);
% 
% % Plot electric field
% figure(2)
% subplot(2,1,1);
% plot(xmesh(x_indices), e_analytic, xmesh(x_indices), real(ebow_freq(idx_plot+np))/dy, xmesh(x_indices), ebow_new(idx_plot+np)/dy);
% xlim([0  3]);
% ylim([-1.5 1.5]);
% legend({'Analytic','Numerical FD', 'Numerical TD'},'Location','southwest');
% title('E polarized in y-direction: E in y-direction')
% 
% % Plot magnetic field
% subplot(2,1,2);
% plot(xmesh(x_indices)-dx/2,h_analytic,xmesh(x_indices)-dx/2,real(hbow_freq(idx_plot+2*np)),xmesh(x_indices)-dx/2,hbow_old(idx_plot+2*np));
% xlim([0  3]);
% ylim([-1.5/Z  1.5/Z]);
% legend({'Analytic','Numerical FD', 'Numerical TD'},'Location','southwest');
% title('E polarized in y-direction: H in z-direction')
% drawnow
% 
% 
% %% Calculate errors for wave with polarization of E in z-direction 
% 
% % Errors electric field
% error_E_FD = norm(real(ebow_freq(idx_plot+np))/dy - e_analytic') / norm(e_analytic);
% error_E_TD = norm(ebow_new(idx_plot+np)/dy - e_analytic') / norm(e_analytic);
% 
% % Errors magnetic field
% error_H_FD = norm(real(hbow_freq(idx_plot+2*np)) - h_analytic') / norm(h_analytic);
% error_H_TD = norm(hbow_old(idx_plot+2*np) - h_analytic') / norm(h_analytic);
% 
% % Display errors:
% disp('Errors for electromagnetic field with polarization of E in y-direction:');
% disp(['Error electric field FD: ', num2str(error_E_FD)]);
% disp(['Error electric field TD: ', num2str(error_E_TD)]);
% disp(['Error magnetic field FD: ', num2str(error_H_FD)]);
% disp(['Error magnetic field TD: ', num2str(error_H_TD)]);
