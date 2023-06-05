% Einfacher Solver für elektrostatische Probleme.
%
%   Eingabe
%   msh         Kanonisches kartesisches Gitter
%   mu          Permeabilität entsprechend Methode createMmu
%   jbow        Stromgitterfluss entsprechend Methode calcHi
%
%   Rückgabe
%   hbow        Magnetische Gitterspannung
%   bbow        Magnetische Gitterfluss
%   relRes      Relatives Residuum für jeden Iterationsschritt

function [hbow, bbow, relRes] = solveMS(msh, mu, jbow)

    bc = [ 0 0 0 0 0 0 ];

    % Erzeugung topologische Matrizen
    [~, ~, st] = createTopMats(msh);

    % Erzeugung geometrische Matrizen
    [ds, ~, da, dat] = createGeoMats(msh);

    % Erzeugung Materialmatrix Permeabilität
    np = msh.np;
    Mmu = createMeps(msh, ds, da, dat, mu, bc);

    % Berechnung Systemmatrix
    A = st * Mmu * st';

    % Berechnen des Hilfsfelds mit calcHi
    hibow = calcHi(msh,jbow);

    % Ladungsvektor berechnen
    qm = st * Mmu * hibow;

    % Gleichungssystem lösen (bitte das Minus vor der magnetischen Ladung beachten)
    [phi, ~, ~, ~, resVec] = gmres(A, -qm, 20, 1e-10, msh.np);
    % Wenn gmres(20) nicht konvergieren würde, probieren Sie bitte bicgstab
    % [phi, ~, ~, ~, resVec] = bicgstab(A, -qm, 1e-10, msh.np);

    relRes = resVec./norm(qm);  %gmres gibt relRes nur für die letzte
                                %Iteration zurueck, so kann man sie jedoch
                                %für alle berechnen.


    %Magnetische Gitterspannung und Gitterfluss berechnen
    hhbow = st' * phi;
    hbow = hhbow + hibow;
    bbow = Mmu * hbow;

end