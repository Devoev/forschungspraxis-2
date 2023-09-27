
function [zedges,n_zedges, deltas] = initMur_2D(msh, open_bc)
%% Description
%
% This function creates the index lists for the boundary z-edges for a 2D
% domain
%
% Input
% msh               2D mesh object
% open_bc           list of four bool vals (if TRUE-> open bc)
%
% Output
% zedges            struct with the indices of the boundary z-edges
% n_zedges          struct with indices of neigbouring z-edges
% deltas            space discretization of the mesh in this direction


%% Function definition

    % Extract mesh parameters
    np = msh.np;
    nx = msh.nx;
    ny = msh.ny;
    Mx = msh.Mx;
    My = msh.My;

    
    % Get selected boundaries
    b1_bool = open_bc(1);
    b2_bool = open_bc(2);
    b3_bool = open_bc(3);
    b4_bool = open_bc(4);

    % Get increments of equidistant mesh
    dx = abs(msh.xmesh(2)-msh.xmesh(1));
    dy = abs(msh.ymesh(2)-msh.ymesh(1));


    % Calculate indices of boundary z edges
    indy = repmat(1:ny,1,1);
    indx = repmat(1:nx,1,1);

    if b1_bool
        % calculate indices for boundary on the right (x = x_max)
        b1_zedges = 1+(nx-1)*Mx+(indy-1)*My + 2*np;
    end
    if b2_bool
        % calculate indices for boundary on the bottom (y = y_min)
        b2_zedges = 1+(indx-1)*Mx + 2*np;
    end
    if b3_bool
        % calculate indices for boundary on the left (x = x_min)
        b3_zedges = 1+(indy-1)*My + 2*np;
    end
    if b4_bool
        % calculate indices for boundary on the top (y = y_max)
        b4_zedges = 1+(indx-1)*Mx+(ny-1)*My + 2*np;
    end


    % Determine neighbouring z-edges
    b1_n_zedges = b1_zedges-msh.Mx;
    b2_n_zedges = b2_zedges+msh.My;
    b3_n_zedges = b3_zedges+msh.Mx;
    b4_n_zedges = b4_zedges-msh.My;

    
    % Save obtained data in struct
    zedges = struct('b1',b1_zedges,'b2', b2_zedges, 'b3', b3_zedges, 'b4', b4_zedges);
    n_zedges = struct('b1', b1_n_zedges, 'b2', b2_n_zedges, 'b3', b3_n_zedges, 'b4', b4_n_zedges);
    deltas = [dx, dy, dx, dy];

end
