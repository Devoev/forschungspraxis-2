function W = projector(dim, ind)
    % PROJECTOR - Create a projection matrix to project onto a subspace
    
    % Input:
    %   dim - Dimension of the space
    %   ind - Indices of the subspace to be removed
    
    % Output:
    %   W - Projection matrix
    
    % Create an identity matrix of size 'dim'
    P = speye(dim);
    
    % Set the elements at indices specified by 'ind' to 0
    P(ind, ind) = 0;
    
    % Create a temporary index array for the complementary subspace
    ind_tmp = 1:dim;
    
    % Remove the indices specified by 'ind' from the temporary index array
    ind_tmp = ind_tmp(~ismember(ind_tmp, ind));
    
    % Extract columns corresponding to the complementary subspace
    Q = P(:, ind_tmp);
    
    % Compute the projection matrix by multiplying P and Q
    W = sparse(P * Q);
end
