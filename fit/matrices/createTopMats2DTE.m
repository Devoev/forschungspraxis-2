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
% CREATE_TOP_MATS_2DTE Creates the topological matrices in the 2D TE case.
% Inputs:
%   msh  - Mesh struct.
% Outputs:
%   c    - Curl matrix.
%   g    - Grad matrix.
%   st   - Dual source matrix.
   
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
% CREATE_PX Creates the Px matrix.
% Inputs:
%   np  - Number of points in the grid.
% Outputs:
%   px  - Px matrix
    
    row=[1:np,1:(np-1)];
    column=[1:np,2:np];
    values=[-ones(1,np),ones(1,np-1)];
    px=sparse(row,column,values,np,np);

end

function py=createPy(np, nx)
% CREATE_PY Creates the Py matrix.
% Inputs:
%   np  - Number of points in the grid.
% Outputs:
%   py  - Py matrix
    
    row=[1:np,1:(np-nx)];
    column=[1:np,1+nx:np];
    values=[-ones(1,np),ones(1,np-nx)];
    py=sparse(row,column,values,np,np);

end
