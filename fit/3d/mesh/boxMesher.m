% Diskretisiert eine als Boxen vorgegebene Verteilung (z.B.
% Materialverteilung oder Potentiale).
%
%   Eingabe
%   msh             Kanonisches kartesisches Gitter
%   boxes           Boxen definiert durch Gitterindizes und Wert
%   defaultvalue    Standardwert für nicht durch Boxen definierte Bereiche
%
%   Rückgabe
%   vec             Vektor der diskretisierten Verteilung (np-by-1)

function [vec] = boxMesher(msh, boxes, defaultvalue)

    if nargin < 3
        defaultvalue = NaN;
    end

    np = msh.np;
    My = msh.My;
    Mz = msh.Mz;
    vec = defaultvalue * ones(np, 1);
    
    for box=boxes
        xmin = box.box(1);
        xmax = box.box(2);
        ymin = box.box(3);
        ymax = box.box(4);
        zmin = box.box(5);
        zmax = box.box(6);
        
        diffx = xmax-xmin+1;
        diffy = ymax-ymin+1;
        diffz = zmax-zmin+1;
        
        inds = xmin:xmax;
        inds = repmat(inds, [diffy 1])' + repmat(((ymin:ymax)-1)* My, [diffx 1]);
        inds = reshape(inds, [1 diffx*diffy]);
        inds = repmat(inds, [diffz 1])' + repmat(((zmin:zmax)-1)* Mz, [diffx*diffy 1]);
        
        vec(inds) = box.value;
    end
    
end