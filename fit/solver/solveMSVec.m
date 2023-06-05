% Aufgabe 2

% Einfacher Solver fuer magnetostatische Probleme der einen Vektor-
% potentialansatz verwendet
%
%   Eingabe
%   msh         Kanonisches kartesisches Gitter
%   mui         Permeabilitaet entsprechend Methode createMmui
%   jbow        Stromgitterfluss entsprechend Methode calcHi
%
%   Rueckgabe
%   hbow        Magnetische Gitterspannung
%   bbow        Magnetische Gitterfluss
%   relRes      Relatives Residuum fuer jeden Iterationsschritt

function [abow, hbow, bbow, relRes] = solveMSVec(msh, mui, jbow)

    % Erzeugung topologische Matrizen
    [c, ~, ~] = createTopMats(msh);

    % Erzeugung geometrische Matrizen
    [ds, dst, da, ~] = createGeoMats(msh);

    % Erzeugung Materialmatrix inverse Permeabilitaet
    mmui = createMmui(msh, ds, dst, da, mui);

    % Berechnung Systemmatrix
    A = c'*mmui*c;

    % Gleichungssystem loesen
    [abow, ~, ~, ~, resVec] = gmres(A, jbow, 20, 1e-8, 500);
    % Wenn gmres(20) nicht konvergieren würde, probieren Sie bitte bicgstab
    % [abow, ~, ~, ~, resVec] = bicgstab(A, jbow, 1e-8, 500);
    relRes = resVec./norm(jbow);    % gmres gibt relRes nur fuer die letzte
                                    % Iteration zurueck, so kann man sie
                                    % jedoch fuer alle berechnen.

    % Magnetische Gitterspannung und Gitterfluss berechnen
    bbow = c*abow;
    hbow = mmui*bbow;

end