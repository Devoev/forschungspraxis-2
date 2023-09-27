%% Gitter erstellen
xmesh = linspace(0,1,61);
ymesh = linspace(0,1,61);

% set open boundary for each side (overrides the initial bc)
open_bc = [true, true, true, true];  % [L1, L2, L3, L4];


msh = cartMesh(xmesh, ymesh); 

Mx = msh.Mx;
My = msh.My;
nx = msh.nx;
ny = msh.ny;
np = msh.np;

%% Erzeugung der topologischen und geometrischen Matrizen
[c, s, st] = createTopMats(msh);
[ds, dst, da, dat] = createGeoMats(msh);

% Erzeugung der Materialmatrizen: Mmui, Mmu, Meps, Mepsi
epsilon = 1;
mui = 1;

Meps = createMeps(msh, ds, dat, epsilon);
Mepsi = nullInv(Meps);

Mmui = createMmui(msh, dst, da, mui);
Mmu = nullInv(Mmui);

% init mur edges
[mur_edges, mur_n_edges, mur_deltas] = initMur(msh, open_bc);

%% Experimentelle Bestimmung mithilfe der Energie des Systems

% Parameter der Zeitsimulation
sigma = 6e-10;
dt = 1e-11;
tend = 10*sigma;
steps = floor(tend/dt);
sourcetype= 1;  % 1: Gauss Anregung, 2: Harmonisch, 3: Konstante Anregung

% Anregung jsbow als anonyme Funktion, die einen 3*np Vektor zurückgibt
% Anregung wird später in Schleife für t>2*sigma_gauss abgeschnitten, also null gesetzt
jsbow_space = zeros(3*np, 1);
x_L = ceil(msh.nx/2);
y_L = ceil(msh.ny/2);


n = 1 + x_L*Mx + y_L*My + 2*np;
jsbow_space(n) = 1;  


jmax = 1;

% Gauss Anregung
jsbow_gauss = @(t)(jsbow_space * jmax * exp(-4*((t-sigma)/sigma)^2));

% Harmonische Anregung (optional)
f = 1e9;
jsbow_harm = @(t)(jsbow_space * jmax * sin(2*pi*f*t));

% Konstante Anregung (optional, für t>0)
jsbow_const = @(t)(jsbow_space * jmax);

% Initialisierungen
ebow_new = sparse(3*np,1);
hbow_new = sparse(3*np,1);
energy = zeros(1,steps);
leistungQuelle = zeros(1,steps);

% Plot parameter für "movie"
figure(1)
%zlimit = 30/(delta_x(1));
zlimit = 700;
draw_only_every = 4;

% Zeitintegration
for ii = 1:steps
    % Zeitpunkt berechnen
    t = ii*dt;

    % alte Werte speichern
    ebow_old = ebow_new;
    hbow_old = hbow_new;

    if sourcetype == 1
        % Stromanregung js aus jsbow(t) für diesen Zeitschritt berechnen
        if t <= 2*sigma
            js = jsbow_gauss(t);
        else
            %js = 0;
            js = sparse(3*np,1);
        end
    elseif sourcetype == 2
        % Harmonische Anregung
        js = jsbow_harm(t);
    elseif sourcetype == 3
        % konstante Anregung
        js = jsbow_const(t);
    end
    
    % Leapfrogschritt durchführen
    [hbow_new,ebow_new] = leapfrog(hbow_old, ebow_old, js, Mmui, Mepsi, c, dt);
    % apply open boundary w. mur cond
    ebow_new = applyMur(mur_edges, mur_n_edges, mur_deltas, ebow_old, ebow_new, dt);
    
    % Feld anzeigen
    if mod(ii, draw_only_every)
        z_plane = 1;
        idx2plot = 2*np+1:3*np;
		ebow_mat = reshape(ebow_new(idx2plot),nx,ny);
    	figure(1)
        mesh(ebow_mat)
        xlabel('i')
        ylabel('j')
        zlabel(['z-Komponente des E-Feldes für z=',num2str(z_plane)])
		axis([1 nx 1 ny -zlimit zlimit])
		caxis([-zlimit zlimit])
        drawnow
    end

    % Gesamtenergie und Quellenenergie für diesen Zeitschritt berechnen
    energy_t = 0.5 * (ebow_new' * Meps*ebow_new + hbow_new' * Mmu * hbow_new);
    leistungQuelle_t = ebow_new' * js;

    % Energiewerte speichern
    energy(ii) =  energy_t;
    leistungQuelle(ii) = leistungQuelle_t;
end

% Anregungsstrom über der Zeit plotten
figure(2)
jsbow_plot = zeros(1,steps);
for step = 1:steps
    if sourcetype == 1
        jsbow_spatial = jsbow_gauss(step*dt);
    elseif sourcetype == 2
        jsbow_spatial = jsbow_harm(step*dt);
    elseif sourcetype == 3
        jsbow_spatial = jsbow_const(step*dt);
    end
    nonzero_idx = find(jsbow_spatial~=0);
    jsbow_plot(step) = jsbow_spatial(nonzero_idx);
end
plot(dt:dt:dt*steps, jsbow_plot);
xlabel('t in s');
ylabel('Anregungsstrom J in A');

% Energie über der Zeit plotten
figure(3); clf;
plot (dt:dt:dt*steps, energy)
legend(['Zeitschritt: ', num2str(dt)])
xlabel('t in s')
ylabel('Energie des EM-Feldes W in J')

% Zeitliche Änderung der Energie (Leistung)
leistungSystem = diff(energy) / dt;
figure(4); clf;
hold on
plot(2*dt:dt:dt*(steps), leistungSystem)
plot(dt:dt:dt*steps, leistungQuelle, 'r')
hold off
legend('Leistung System', 'Leistung Quelle')
xlabel('t in s')
ylabel('Leistung P in W')