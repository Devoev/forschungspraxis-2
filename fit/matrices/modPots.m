%   Modifiziert die Systemmatrix und den Ladungsvektor mit
%   vorgegebenen Potentialen.
%
%   Eingabe
%   A               Systemmatrix
%   q               Ladungsvektor
%   pots            Potentialvektor
%
%   Rückgabe
%   A               Modifizierte Systemmatrix
%   q               Modifizierter Ladungsvektor

function [A, q] = modPots(A, q, pots)

    pots = reshape(pots, [numel(pots) 1]);
    q = reshape(q, [numel(q) 1]);
    eqsinds = isnan(pots);
    bcinds = ~eqsinds;
    q = q(eqsinds) - A(eqsinds, bcinds) * pots(bcinds);
    A = A(eqsinds,eqsinds);
    
end