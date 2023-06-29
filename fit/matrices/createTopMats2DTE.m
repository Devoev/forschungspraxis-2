% Erzeugt die topologischen Matrizen für ein kanonisches Gitter.
%
%   Eingabe
%   msh         Kanonisches Gitter
%
%   Rückgabe
%   c           Curl-Matrix
%   s           Source-Matrix
%   st          Duale Source-Matrix

function [c, g, st]=createTopMats2DTE(msh)
   
    np=msh.np;
    px=createPx(np);
    py=createPy(np,msh.nx);
    
    % Curl-Matrix
    c = [py; -px];

    % Grad Matrix
    g = [px; py];
    
    % Source-Matrix
    st = [-px', -py'];

end

function px=createPx(np)
    
    row=[1:np,1:(np-1)];
    column=[1:np,2:np];
    values=[-ones(1,np),ones(1,np-1)];
    px=sparse(row,column,values,np,np);

end

function py=createPy(np, nx)
    
    row=[1:np,1:(np-nx)];
    column=[1:np,1+nx:np];
    values=[-ones(1,np),ones(1,np-nx)];
    py=sparse(row,column,values,np,np);

end
