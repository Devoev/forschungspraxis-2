
function vec = applyMur_2D(zedges, n_zedges, deltas, vec_old, vec_new, dt)
%% Description
%
% Creates a basic 2D mesh object
%
% Input
% zedges            outer z edges for open boundary
% n_zedges          inner z edges for open boundary 
% deltas            increments in each direction
% vec_old           old vector of primary z component
% vec_new           new vector of primary z component
% dt                time step size
%
% Output
% vec               new result vector with effective open boundary


%% Function definition

    % Initialize new vector
    vec = sparse(vec_new);
    
    % boundary1
    edge = zedges.b1;
    n_edge = n_zedges.b1;
    delta = deltas(1);
    vec(edge) = mur(edge, n_edge, delta, dt, vec_old, vec_new);
    
    % boundary2
    edge = zedges.b2;
    n_edge = n_zedges.b2;
    delta = deltas(2);
    vec(edge) = mur(edge, n_edge, delta, dt, vec_old, vec_new);
    
    % boundary3
    edge = zedges.b3;
    n_edge = n_zedges.b3;
    delta = deltas(3);
    vec(edge) = mur(edge, n_edge, delta, dt, vec_old, vec_new);
    
    % boundary4
    edge = zedges.b4;
    n_edge = n_zedges.b4;
    delta = deltas(4);
    vec(edge) = mur(edge, n_edge, delta, dt, vec_old, vec_new);


    function edgevals = mur(edge, n_edge, delta, dt, vec_old, vec_new)
        c0 = 3e8;
        edgevals = vec_old(n_edge) + ((c0*dt-delta)/(c0*dt+delta)).*(vec_new(n_edge)-vec_old(edge));
    end

end
