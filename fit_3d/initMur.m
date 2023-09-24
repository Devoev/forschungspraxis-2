function [zedges,n_zedges, deltas] = initMur(msh, open_bc)
%CREATE_MUR This function creates the index lists for the boundary edges
%   input:  msh - geometry mesh of the domain
%           open_bc - list of four bool vals (if TRUE-> open bc)
%
%   output: zedges - struct with the indices of the boundary z-edges
%           n_zedges - struct with indices of neigbouring z-edges
%           deltas - space discretization of the mesh in this direction
%

b1_bool = open_bc(1);
b2_bool = open_bc(2);
b3_bool = open_bc(3);
b4_bool = open_bc(4);

dx = abs(msh.xmesh(2)-msh.xmesh(1));
dy = abs(msh.ymesh(2)-msh.ymesh(1));

indices = 1:msh.Mz;
% get element indices of selected boundaries b1-b4
b1 = [];
b2 = [];
b3 = [];
b4 = [];

if b1_bool
    b1 = indices(~(mod(indices,msh.nx)));   % y right edges NIX??
end
if b2_bool
    b2 = 1:msh.nx;                          % X bottom edges
end
if b3_bool
    b3 = indices(~(mod(indices,msh.nx)-1)); % y left edges
end
if b4_bool
    b4 = (msh.Mz-msh.nx+1):msh.Mz;          % top X edges
end

% build z-edge indices for all z-layers
b1_edges = [];
b2_edges = [];
b3_edges = [];
b4_edges = [];
for z_idx = 1:msh.nz
    b1_edges = [b1_edges, b1+(z_idx-1)*msh.Mz];
    b2_edges = [b2_edges, b2+(z_idx-1)*msh.Mz];
    b3_edges = [b3_edges, b3+(z_idx-1)*msh.Mz];
    b4_edges = [b4_edges, b4+(z_idx-1)*msh.Mz];
end

% add z-edge offset
b1_zedges = b1_edges+2*msh.np;
b2_zedges = b2_edges+2*msh.np;
b3_zedges = b3_edges+2*msh.np;
b4_zedges = b4_edges+2*msh.np;

% get neighbouring z-edges
b1_n_zedges = b1_zedges-msh.Mx;
b2_n_zedges = b2_zedges+msh.My;
b3_n_zedges = b3_zedges+msh.Mx;
b4_n_zedges = b4_zedges-msh.My;


zedges = struct('b1',b1_zedges,'b2', b2_zedges, 'b3', b3_zedges, 'b4', b4_zedges);
n_zedges = struct('b1', b1_n_zedges, 'b2', b2_n_zedges, 'b3', b3_n_zedges, 'b4', b4_n_zedges);
deltas = [dx, dy, dx, dy];
end

