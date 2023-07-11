%Erzeugt die Pseudo-Inverse einer diagonalen Matrix A, d.h. invertiert
%jeden Diagonaleintrag ungleich null.
%
%   Eingabe
%   A           Zu invertierende Matrix
%
%   Rückgabe
%   AInv        Pseudo-Inverse der Matrix

function AInv = nullInv(A)

    [indm, indn, values] = find(A);
    AInv = sparse(indm, indn, 1./values, size(A,1), size(A,2));

end

