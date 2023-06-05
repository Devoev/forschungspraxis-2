%   Einfacher Solver für elektrostatische Probleme
%
%   Eingabe
%   msh         Kanonisches kartesisches Gitter
%   epsilon     Permittivität entsprechend Methode createMeps
%   pots        Potentiale entsprechend Methode boxMesher
%   q           Ladungsvektor entsprechend Methode boxMesher
%
%   Rückgabe
%   phi         Potentialvektor Lösung
%   ebow        Elektrische Gitterspannung
%   dbow        Elektrischer Gitterfluss
%   relRes      Relatives Residuum für jeden Iterationsschritt

function [phi, ebow, dbow, relRes] = solveES(msh, epsilon, pots, q)

    % Erzeugung topologische Matrizen
    [~, ~, st] = createTopMats(msh);

    % Erzeugung geometrische Matrizen
    [ds, ~, da, dat] = createGeoMats(msh);

    % Erzeugung der Materialmatrix der Permittivität mit createMeps
    Meps = createMeps(msh, ds, da, dat, epsilon);

    % Berechnung Systemmatrix
    A = st*Meps*st';

    % Modifikation Systemmatrix und Ladungsvektor mit modPots
    [A, q] = modPots(A, q, pots);

    % Gleichungssystem lösen mit gmres(20)
    [x, ~, ~, ~, resVec] = gmres(A, q, 20, 1e-13, msh.np);
    % Wenn gmres(20) nicht konvergieren würde, probieren Sie bitte bicgstab
    % [x, ~, ~, ~, resVec] = bicgstab(A, q, 1e-13, msh.np);

    % Bestimmung des Residuums
    relRes = resVec./norm(q);   % gmres gibt relRes nur für die letzte
                                % Iteration zurück, so kann man sie jedoch
                                % für alle berechnen.

    % phi aus x bestimmen (eingeprägte Potentiale wieder einfügen)
    idx = isnan(pots);
    phi = pots;
    phi(idx) = x;

    % Elektrische Gitterspannung und Gitterfluss berechnen
    ebow = st'*phi;
    dbow = Meps*ebow;

end