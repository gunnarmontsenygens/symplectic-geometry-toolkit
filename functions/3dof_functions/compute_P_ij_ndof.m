function P_ij_mtx = compute_P_ij_ndof(n, i, j)
%==========================================================================
%
% Constructs the permutation matrix P_{ij} \in \mathbb{R}^{n \times n}
% that swaps entries i and j of a vector through left multiplication.
%
% MATRIX DESCRIPTION:
% The matrix P_{ij} is an identity matrix with rows i and j exchanged.
% For a vector x \in \mathbb{R}^n:
%
%   x_swapped = P_{ij} * x
%
% produces a new vector in which:
%
%   - Entry i of x is moved to position j
%   - Entry j of x is moved to position i
%   - All other entries remain unchanged
%
% PROPERTIES:
%
%   - P_{ij}^{-1} = P_{ij}^{T} = P_{ij}
%   - P_{ij}^2 = I
%   - det(P_{ij}) = -1
%
% Author: G. Montseny
% Date: May 28, 2026
%
% INPUT:               Description                             Units
%
%   n           -   Dimension of the vector space             [-]
%   i           -   First index to swap                       [-]
%   j           -   Second index to swap                      [-]
%
% OUTPUT:              Description                             Units
%
%   P_ij_mtx    -   Permutation matrix exchanging             [-]
%                   entries i and j
%
%==========================================================================

    % Preallocate matrix
    P_ij_mtx = zeros(n, n);

    % Build matrix
    for k = 1 : n
        if k == i
            P_ij_mtx(k, :) = [zeros(1,j-1), 1, zeros(1, n-j)];
        elseif k == j
            P_ij_mtx(k, :) = [zeros(1, i-1), 1, zeros(1, n-i)];
        else
            P_ij_mtx(k, :) = [zeros(1, k-1), 1, zeros(1, n-k)];
        end
    end
end