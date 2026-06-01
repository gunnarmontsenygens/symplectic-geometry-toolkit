function Psi_tilde_mtx = compute_Psi_tilde_3dof(t0, tf, x0_vec, xf_vec, Phi_mtx, params)
%==========================================================================
%
% Computes the transformed State Transition Matrix (STM) in the
% coordinate-momentum pair ordering:
%
%   Psi_tilde = P_x2z * Psi * P_x2z^T
%
% for a spatial 3-DoF Hamiltonian dynamical system using the
% ortho-symplectic transformation introduced by Scheeres & Boodram (2025).
%
% MODEL DESCRIPTION:
% The transformed STM Psi is first constructed using the orthogonal
% symplectic bases evaluated at the initial and final states:
%
%   Psi = R_f^T * Phi * R_o
%
% where:
%
%   Phi    - State Transition Matrix (STM)
%   R_o    - Initial ortho-symplectic transformation matrix
%   R_f    - Final ortho-symplectic transformation matrix
%
% The matrix Psi is naturally expressed in the canonical ordering:
%
%   [q_1 q_2 q_3 p_1 p_2 p_3]
%
% This function then applies a permutation transformation to obtain:
%
%   Psi_tilde
%
% in the coordinate-momentum pair ordering:
%
%   [q_1 p_1 q_2 p_2 q_3 p_3]
%
% which exposes the block structure described in the reference paper.
%
% The permutation is performed using:
%
%   z = P_x2z x
%
% and:
%
%   Psi_tilde = P_x2z Psi P_x2z^T
%
% CURRENT IMPLEMENTATION:
%   - Spatial 3-DoF systems only
%   - Full state vectors and STM are provided directly
%   - Supports both Lagrangian and Hamiltonian formulations
%
% STATE DEFINITION:
% Full state vectors must be provided as:
%
%   Lagrangian:
%
%       x_vec = [x y z vx vy vz]^T
%
%   Hamiltonian:
%
%       x_vec = [x y z px py pz]^T
%
% STM DEFINITION:
% The STM must correspond to the same formulation as the state vector:
%
%       Phi = d x_f / d x_0
%
% FORMULATION HANDLING:
% If the STM is provided in Lagrangian coordinates, it is internally
% transformed into canonical Hamiltonian coordinates before constructing
% Psi and Psi_tilde.
%
% Author: G. Montseny
% Date: May 29, 2026
%
% INPUTS:                   Description                          Units
%
%   t0              - Initial time                              [-]
%   tf              - Final time                                [-]
%   x0_vec          - Initial state vector (6x1)                [-]
%   xf_vec          - Final state vector (6x1)                  [-]
%   Phi_mtx         - State Transition Matrix (6x6)             [-]
%   params          - Parameter struct                          [-]
%
% OUTPUTS:                  Description                          Units
%
%   Psi_tilde_mtx   - Transformed STM in paired ordering        [-]
%
%==========================================================================

    % ---------------------------------------------------------------------
    %              COMPUTE PSI
    % ---------------------------------------------------------------------

    Psi_mtx = compute_Psi_3dof(t0, tf, x0_vec, xf_vec, Phi_mtx, params);

    % ---------------------------------------------------------------------
    %              COMPUTE P_x2z
    % ---------------------------------------------------------------------

    P_x2z_mtx = compute_P_x2z_3dof();

    % ---------------------------------------------------------------------
    %              COMPUTE PSI_TILDE
    % ---------------------------------------------------------------------

    Psi_tilde_mtx = P_x2z_mtx * Psi_mtx * P_x2z_mtx';

end