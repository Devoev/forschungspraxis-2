%   Erstellt einen Höhenlinienplot eines Potentials
%
%   Eingabe
%   msh             Kanonisches kartesisches Gitter
%   phi             Potentialvektor
%   indz            z-Index der darzustellenden x-y-Schnittebene
%
%   Rückgabe
%   figure(1)       Plot des Potentials (wird abgespeichert in plotPot.pdf)

function plotPotential(msh, phi, indz)

    pick = ((indz-1)*msh.nx*msh.ny+1):indz*msh.nx*msh.ny;
    [X, Y] = meshgrid(msh.xmesh, msh.ymesh);
    figure(1)
        contourf(X, Y, reshape(phi(pick), msh.nx, msh.ny)', 15);
        xlabel('x')
        ylabel('y')
end