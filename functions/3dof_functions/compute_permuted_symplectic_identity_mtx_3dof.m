function J_p_mtx = compute_permuted_symplectic_identity_mtx_3dof()
%==========================================================================
%
% Computes the permuted symplectic identity matrix for a 3-DoF Hamiltonian
% system in coordinate-momentum pair ordering.
%
% MODEL DESCRIPTION:
% The canonical symplectic identity matrix is naturally defined in the
% standard Hamiltonian state ordering:
%
%   x = [q_1 q_2 q_3 p_1 p_2 p_3]^T
%
% where the symplectic structure is represented by:
%
%   J = [ 0   I ]
%       [ -I  0 ]
%
% In some applications, it is convenient to work in the coordinate-momentum
% pair ordering:
%
%   z = [q_1 p_1 q_2 p_2 q_3 p_3]^T
%
% The corresponding symplectic identity matrix in the permuted basis is:
%
%   J_p = P_x2z * J * P_x2z^T
%
% where P_x2z is the permutation matrix that maps the canonical ordering
% into the coordinate-momentum pair ordering.
%
% This matrix preserves the symplectic structure in the permuted
% coordinates and can be used to evaluate symplecticity conditions such as:
%
%   Phi^T * J_p * Phi = J_p
%
% for state transition matrices expressed in the permuted basis.
%
% Author: G. Montseny
% Date: May 29, 2026
%
% OUTPUT:              Description                                   Units
%
%  J_p_mtx    -   permuted symplectic identity matrix                [-]
%
%==========================================================================

    % Compute permutation matrices
    P_x2z_mtx = compute_P_x2z_3dof();
    P_z2x_mtx = P_x2z_mtx';

    % Compute non-permuted symplectic identity matrix
    J_mtx = compute_symplectic_identity_mtx(3);

    % Compute permuted symplectic identity matrix
    J_p_mtx = P_x2z_mtx * J_mtx * P_z2x_mtx;

end