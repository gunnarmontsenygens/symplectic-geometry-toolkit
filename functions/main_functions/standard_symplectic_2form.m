function w = standard_symplectic_2form(xi1_vec, xi2_vec)
%==========================================================================
%
% Computes the standard symplectic 2-form between two vectors in R^{2n}:
%
%   w(xi1,xi2) = xi1^T * J * xi2
%
% where J is the canonical symplectic identity matrix.
%
% Author: G. Montseny
% Date: May 20, 2026
%
% INPUTS:                   Description                          Units
%
%   xi1_vec      - First vector in R^{2n}                       [-]
%   xi2_vec      - Second vector in R^{2n}                      [-]
%
% OUTPUTS:                  Description                          Units
%
%   w            - Symplectic 2-form value                      [-]
%
%==========================================================================
        

    % Initialization
    xi1_vec = xi1_vec(:);
    xi2_vec = xi2_vec(:);

    if length(xi1_vec) ~= length(xi2_vec)
        error('Input vectors must have the same dimension.')
    end

    N = length(xi2_vec);
    n = N/2;


    % Compute symplectic identity matrix
    J_mtx = compute_symplectic_identity_mtx(n);

    % Output
    w = xi1_vec' * J_mtx * xi2_vec;

end