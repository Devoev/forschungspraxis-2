function [e_mur] = applyMur_2D(mur_edges, mur_n_edges, mur_deltas, e_old, e_new, dt, bc)
% applyMur_2D applies open boundary condition for time domain for a 2D mesh
%
% Inputs:
%   mur_edges       - Outer edges for open boundary
%   mur_n_edges     - Inner edges for open boundary
%   mur_deltas      - Increments in each direction
%   vec_old         - Old vector of primary component
%   vec_new         - New vector of primary componen
%   dt              - Time step size
%   bc              - Boundary condition struct
%
% Outputs:
%   ebow    - Integrated electric field
%             Unit in x- and y-direction: V; Unit in z-direction: V/m;
%             Size (3*np, 1)
%   hbow    - Integrated magnetic field
%             Unit in x- and y-direction: A; Unit in z-direction: A/m;
%             Size (3*np, 1)


% Calculate speed of light based on constants
if any(ismember(fieldnames(bc),"epsilon0"))
    c0 = sqrt(bc.mu0i/bc.epsilon0);
else
    c0 = 3e8;
end

% Initialize new vector
e_mur = e_new;

%% Z-direction ebow border edges
% boundary1
edge = mur_edges.b1_z;
n_edge = mur_n_edges.b1_z;
delta = mur_deltas(1);
e_mur(edge) = mur(edge, n_edge, delta, dt, e_old, e_new, c0);
% boundary2
edge = mur_edges.b2_z;
n_edge = mur_n_edges.b2_z;
delta = mur_deltas(2);
e_mur(edge) = mur(edge, n_edge, delta, dt, e_old, e_new, c0);
% boundary3
edge = mur_edges.b3_z;
n_edge = mur_n_edges.b3_z;
delta = mur_deltas(3);
e_mur(edge) = mur(edge, n_edge, delta, dt, e_old, e_new, c0);
% boundary4
edge = mur_edges.b4_z;
n_edge = mur_n_edges.b4_z;
delta = mur_deltas(4);
e_mur(edge) = mur(edge, n_edge, delta, dt, e_old, e_new, c0);

%% X-direction ebow border edges
% boundary1
edge = mur_edges.b2_x;
n_edge = mur_n_edges.b2_x;
delta = mur_deltas(2);
e_mur(edge) = mur(edge, n_edge, delta, dt, e_old, e_new, c0);
% boundary4
edge = mur_edges.b4_x;
n_edge = mur_n_edges.b4_x;
delta = mur_deltas(2);
e_mur(edge) = mur(edge, n_edge, delta, dt, e_old, e_new, c0);
%% Y-direction ebow border edges
% boundary1
edge = mur_edges.b1_y;
n_edge = mur_n_edges.b1_y;
delta = mur_deltas(1);
e_mur(edge) = mur(edge, n_edge, delta, dt, e_old, e_new, c0);
% boundary4
edge = mur_edges.b3_y;
n_edge = mur_n_edges.b3_y;
delta = mur_deltas(1);
e_mur(edge) = mur(edge, n_edge, delta, dt, e_old, e_new, c0);



function edgevals = mur(edge, n_edge, delta, dt, vec_old, vec_new, c0)
    edgevals = vec_old(n_edge) + ((c0*dt-delta)/(c0*dt+delta)).*(vec_new(n_edge)-vec_old(edge));
end


end
