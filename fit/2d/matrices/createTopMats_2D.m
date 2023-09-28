
function [c, s, st] = createTopMats_2D(msh)
%% Description
%
% Create curl and source matrices for 2D case
%
% Input
% msh               2D mesh object
%
% Output
% c                 primary curls matrix
% s                 primary source matrix
% st                dual source matrix


%% Function definition
   
    np=msh.np;
    px=createPx(np);
    py=createPy(np,msh.nx);
    
    % Curl-Matrix
    zeroMat=sparse(np,np);
    c=[zeroMat, zeroMat, py; zeroMat, zeroMat, -px; -py, px, zeroMat];
    
    % Source-Matrix
    s=[px, py, zeroMat];
    
    % Duale Source-Matrix
    st=[-px', -py', zeroMat];

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