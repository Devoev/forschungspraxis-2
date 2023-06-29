% Erzeugt die topologischen Matrizen für ein kanonisches Gitter.
%
%   Eingabe
%   msh         Kanonisches Gitter
%
%   Rückgabe
%   c           Curl-Matrix
%   s           Source-Matrix
%   st          Duale Source-Matrix

function [c, s, st]=createTopMats2D(msh)
   
    np=msh.np;
    px=createPx(np);
    py=createPy(np,msh.nx);
    
    % Curl-Matrix
    zeroMat=sparse(np,np);
    c=[zeroMat, -pz, py; pz, zeroMat, -px; -py, px, zeroMat];
    
    % Source-Matrix
    s=[px, py, pz];
    
    % Duale Source-Matrix
    st=[-px', -py', -pz'];

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
