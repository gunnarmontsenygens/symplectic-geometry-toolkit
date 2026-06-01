function J_mtx = compute_symplectic_identity_mtx(n)
%==========================================================================
%
% Computes the canonical symplectic identity matrix:
%
%              [  0   I ]
%   J_mtx  =   [       ]
%              [ -I   0 ]
%
% for a 2n-dimensional Hamiltonian system.
%
% Author: G. Montseny
% Date: May 20, 2026
%
% INPUTS:                   Description                          Units
%
%   n             - Dimension of the configuration space        [-]
%
% OUTPUTS:                  Description                          Units
%
%   J_mtx         - Canonical symplectic identity matrix        [-]
%                   of size [2n x 2n]
%
%==========================================================================

    J_mtx = [zeros(n), eye(n);
            -eye(n), zeros(n)];
end