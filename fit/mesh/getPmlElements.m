function idxPmlElements = getPmlElements(msh, n_layers)
% getPmlElements(msh, n_layers) returns the canonical indices (2D) of all 
% elements which are a part of the PML boundary w. layer count n_layers:
% 
%     + + + + + +
%     + + + + + +     n_layers = 2
%     + + . . + +  
%     + + . . + +     . regular elements
%     + + + + + +     + PML elements
%     + + + + + +
%   
% Input:
%   msh             msh struct of the geometry
%   n_layers        number of PML-Element layers
%
% Output:
%   idxPmlElements   column vector that contains all PML elements' indices
%

n_elemx = msh.nx-1;
n_elemy = msh.ny-1;
n_elems = n_elemx*n_elemy;
elems = (1:1:n_elems)';

% check if requested layers are possible:
if (n_elemy <= n_layers*2) || (n_elemx <= n_layers*2) || n_layers == 0
    f = errordlg("illegal choice of PML layers");
    idxPmlElements = 0;
else
    % get all indices of bottom row(s)
    idx_plm = elems(1:n_elemx*n_layers);
    
    for indy = n_layers:1:n_elemy-n_layers-1
        % append indices of left column(s) for row indy
        idx_plm = [idx_plm; elems(n_elemx*indy+1:n_elemx*indy+n_layers)];
        % append indices of right column(s) for row indy
        idx_plm = [idx_plm; elems(n_elemx*(indy+1)-n_layers+1:n_elemx*(indy+1))];
    end
    % append all indices of top row(s)
    idxPmlElements = [idx_plm; elems(n_elems-n_elemx*n_layers+1:n_elems)];
end