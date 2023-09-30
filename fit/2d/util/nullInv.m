
function AInv = nullInv(A)
%% Description
%
% Function to determine the pseudo inverse to a square diagonal matrix by
% inverting each element not equal to zero
%
% Input
% A                 square diagonal matrix to invert 
%
% Output
% AInv              pseudo inverse of A (square diagonal)


%% Function definition

    % Find indices and elements unequal to zero
    [indm, indn, values] = find(A);

    % Calculate pseudo inverse
    AInv = sparse(indm, indn, 1./values, size(A,1), size(A,2));

end
