% Versuch 6

clear;

%% Gitter erstellen (nicht mesh nennen, da dies ein Matlab-Befehl ist)
nx = 41;
ny = 41;
nz = 2;
xmesh = linspace(0,1,nx);
ymesh = linspace(0,1,ny);
zmesh = linspace(0,1,nz);
msh = cartMesh(xmesh, ymesh, zmesh); 

% Gitterweiten in x-,y- und z-Richtung (äquidistant vorausgesetzt)
delta_x = diff(xmesh);
delta_y = diff(ymesh);
delta_z = diff(zmesh);

Mx = msh.Mx;
My = msh.My;
Mz = msh.Mz;
np = msh.np;
idxSource = 2*np+1+floor(nx/2)*Mx+floor(ny/2)*My;
ind_h = [2*np+1:3*np];
ind_e = [1:2*np,idxSource];
%ind_e = [idxSource];
%ind_h = [];
We = projector(3*np,ind_e);
Wh = projector(3*np,ind_h);
%% Erzeugung der topologischen und geometrischen Matrizen
[c, s, st] = createTopMats(msh);
[ds, dst, da, dat] = createGeoMats(msh);

% Erzeugung der Materialmatrizen: Mmui, Mmu, Meps, Mepsi
% Achtung bei createMmui und createMeps Randbedingungen übergeben, da
% explizites Verfahren und damit Randbedingungen am besten in den
% Materialmatrizen gesetzt werden
eps0 = 8.854e-12;
eps_r = 1;
epsilon = eps0*eps_r;
mu0 = 4e-7*pi;
mu_r = 1;
mu = mu_r*mu0;
mui = 1/mu;
bcs = [1,1,1,1,1,1];

Mmui = createMmui(msh, ds, dst, da, mui, bcs);
Mmu = nullInv(Mmui);

Meps = createMeps(msh, ds, da, dat, epsilon, bcs);
Mepsi = nullInv(Meps);


%% CFL-Bedingung

% Minimale Gitterweite bestimmen
delta_s = min([delta_x delta_y delta_z]);

% Berechnung und Ausgabe des minimalen Zeitschritts mittels CFL-Bedingung
deltaTmaxCFL = sqrt(mu*epsilon*(1/((1/(min(delta_x)^2))+(1/min(delta_y)^2)+(1/(min(delta_z)^2)))));
fprintf('Nach CFL-Bedingung: deltaTmax = %e\n',deltaTmaxCFL);

%% Experimentelle Bestimmung mithilfe der Energie des Systems

% Parameter der Zeitsimulation
sigma = 6e-10;
dt = 5e-11;
tend = 100*sigma;
timeAxis = 0:dt:tend;
Nt = length(timeAxis);
sourcetype= 1;  % 1: Gauss Anregung, 2: Harmonisch, 3: Konstante Anregung

% Anregung jsbow als anonyme Funktion, die einen 3*np Vektor zurückgibt
% Anregung wird später in Schleife für t>2*sigma_gauss abgeschnitten, also null gesetzt

jsbowAmplitude = 1;
jsbow_space = zeros(3*np, 1);
jsbow_space(idxSource) = jsbowAmplitude;

% Gauss Anregung
jsbow_gauss = @(t) (jsbow_space * exp(-4*((t-sigma)/sigma)^2));

% Harmonische Anregung (optional)
jsbow_harm = @(t) (jsbow_space*sin(2*pi*3/tend*t));

% Konstante Anregung (optional, für t>0)
jsbow_const = @(t)(jsbow_space);
% Initialisierungen
edof = length(nonzeros(We))
hdof = length(nonzeros(Wh))
ebow_veryold = full(sparse(edof,1));
ebow_new = full(sparse(edof,1));
hbow_new = full(sparse(hdof,1));
%ebow = full(sparse(edof,Nt));
%hbow = full(sparse(hdof,Nt));
energy = zeros(1,Nt);
leistungQuelle = zeros(1,Nt);


z_plane = 1;
projector_changed = 1;
% Zeitintegration
for ii = 1:Nt
    % Zeitpunkt berechnen
    t = timeAxis(ii);

    % alte Werte speichern
    if t < 4 * sigma
        egiven = jsbow_gauss(t);
    else 
        egiven = zeros(3*np,1)
    end
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    t = timeAxis(ii)
	jsbow = 0*egiven;
    % Leapfrogschritt durchführen
    [hbow_new,ebow_new] = leapfrog(hbow_old,ebow_old,egiven,0*egiven, Mmui,Mepsi,c,dt,We,Wh,jsbow);

	%ebow(:,ii) = ebow_new;
	%hbow(:,ii) = hbow_new;
    
    % Feld anzeigen
    ebowplot = We * ebow_new;
    idx2plot = (2*np+1+(z_plane-1)*Mz):(2*np+z_plane*Mz);
	ebow_mat = reshape(ebowplot(idx2plot),nx,ny);
    figure(1)
    mesh(ebow_mat);
	axis([1 nx 1 ny -1 1])
    drawnow;
    % Zeit zwischen zwei Frames
    
    %sol{ii} = ebow_mat;
    display(ii);
end