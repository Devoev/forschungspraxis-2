
function vec = boxMesher_2D(msh, boxes, defaultvalue)
%% Description
%
% Function to assign material constants to regions on the 2D calculation
% domain defined in boxes
%
% Input
% msh               2D mesh object
% boxes             struct, which contains the information about the
%                   regions:
%                   box = [x_min, x_max, y_min, y_max]; (Indices, not val!)
%                   value = value of material in the region
%                   
%                   Example:    (Based on 4x4 point domain)
%                               boxes(1).box = [1, 2, 1, 4];
%                               boxes(1).value = 2;
%                               boxes(2).box = [2, 4, 1, 4];
%                               boxes(2).value = 1;
%
% Output
% vec               Vector containing the material parameter to each volume
%                   on the mesh


%% Function definition

    np = msh.np;
    nx = msh.nx;
    ny = msh.ny;
    Mx = msh.Mx;
    My = msh.My;
    
    % Set up vector with default value for each cell
    vec = defaultvalue * ones(np, 1);
    
    % Iterate over each material area
    for box = boxes
        
        % Get boundary indices
        xmin = box.box(1);
        xmax = box.box(2);
        ymin = box.box(3);
        ymax = box.box(4);
 
        % Determine indices of all volumes of the region
        inds = (xmin:xmax-1) + (ymin-1) * My;
        inds = repmat(inds, ymax-ymin, 1) + ones(1,xmax-xmin) .* ((ymin:ymax-1) - ymin)' * My;
        inds = reshape(inds', 1, (xmax-xmin)*(ymax-ymin));
        
        % Assign indices
        vec(inds) = box.value;

    end

    % Set values in ghost volumes to zero
    gvx = 1 + (nx-1)*Mx + (0:ny-1) * My;
    gvy = 1 + (0:nx-1)*Mx + (ny-1) * My;
    vec(gvx) = 0;
    vec(gvy) = 0;

end
