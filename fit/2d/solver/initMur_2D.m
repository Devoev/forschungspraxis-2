function [zedges,n_zedges, deltas] = initMur_2D(msh, bc)
% initMur_2D creates the index lists for the boundary z-edges for a 2D
% domain 
%
% Inputs:
%   msh             - Mesh struct
%   bc              - Boundary condition struct
%
% Outputs:
%   zedges          - struct with the indices of the boundary z-edges
%   n_zedges        - struct with indices of neigbouring z-edges
%   deltas          - space discretization of the mesh in this direction


% Extract mesh parameters
np = msh.np;
nx = msh.nx;
ny = msh.ny;
Mx = msh.Mx;
My = msh.My;


% Get selected boundaries
b1_bool = bc.open_bc(1);
b2_bool = bc.open_bc(2);
b3_bool = bc.open_bc(3);
b4_bool = bc.open_bc(4);

% Get increments of equidistant mesh
dx = abs(msh.xmesh(2)-msh.xmesh(1));
dy = abs(msh.ymesh(2)-msh.ymesh(1));


% Calculate indices of boundary z edges
indy = repmat(1:ny,1,1);
indx = repmat(1:nx,1,1);

b1_zedges = [];
b2_zedges = [];
b3_zedges = [];
b4_zedges = [];

b1_yedges = [];
b2_yedges = [];
b3_yedges = [];
b4_yedges = [];

b1_xedges = [];
b2_xedges = [];
b3_xedges = [];
b4_xedges = [];

if b2_bool
    % calculate indices for boundary on the right (x = x_max)
    b1_xedges = 1+(nx-1)*Mx+(indy-1)*My;
    b1_yedges = b1_xedges + np;
    b1_zedges = b1_xedges + 2*np;
end
if b1_bool
    % calculate indices for boundary on the bottom (y = y_min)
    b2_xedges = 1+(indx-1)*Mx;
    b2_yedges = b2_xedges + np;
    b2_zedges = b2_xedges + 2*np;
end
if b4_bool
    % calculate indices for boundary on the left (x = x_min)
    b3_xedges = 1+(indy-1)*My;
    b3_yedges = b3_xedges + np;
    b3_zedges = b3_xedges + 2*np;

end
if b3_bool
    % calculate indices for boundary on the top (y = y_max)
    b4_xedges = 1+(indx-1)*Mx+(ny-1)*My;
    b4_yedges = b4_xedges + np;
    b4_zedges = b4_xedges + 2*np;
end


% Determine neighbouring z-edges
b1_n_zedges = b1_zedges-msh.Mx;
b2_n_zedges = b2_zedges+msh.My;
b3_n_zedges = b3_zedges+msh.Mx;
b4_n_zedges = b4_zedges-msh.My;

b1_n_yedges = b1_n_zedges - np;
b2_n_yedges = b2_n_zedges - np;
b3_n_yedges = b3_n_zedges - np;
b4_n_yedges = b4_n_zedges - np;

b1_n_xedges = b1_n_zedges - 2*np;
b2_n_xedges = b2_n_zedges - 2*np;
b3_n_xedges = b3_n_zedges - 2*np;
b4_n_xedges = b4_n_zedges - 2*np;


% Save obtained data in struct
zedges = struct('b1_z',b1_zedges,'b2_z', b2_zedges, 'b3_z', b3_zedges, 'b4_z', b4_zedges, ...
    'b1_y',b1_yedges,'b2_y', b2_yedges, 'b3_y', b3_yedges, 'b4_y', b4_yedges, ...
    'b1_x',b1_xedges,'b2_x', b2_xedges, 'b3_x', b3_xedges, 'b4_x', b4_xedges);
n_zedges = struct('b1_z', b1_n_zedges, 'b2_z', b2_n_zedges, 'b3_z', b3_n_zedges, 'b4_z', b4_n_zedges, ...
    'b1_y', b1_n_yedges, 'b2_y', b2_n_yedges, 'b3_y', b3_n_yedges, 'b4_y', b4_n_yedges, ...
    'b1_x', b1_n_xedges, 'b2_x', b2_n_xedges, 'b3_x', b3_n_xedges, 'b4_x', b4_n_xedges);
deltas = [dx, dy, dx, dy];


end
