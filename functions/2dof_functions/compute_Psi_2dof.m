function Psi_mtx = compute_Psi_2dof(t0, tf, x0_vec, xf_vec, Phi_mtx, params)

%==========================================================================
%
% Computes the transformed State Transition Matrix (STM):
%
%   Psi = R_f^T * Phi * R_o
%
% for a planar 2-DoF Hamiltonian subsystem using the ortho-symplectic
% transformation introduced by Scheeres & Boodram (2025).
%
% MODEL DESCRIPTION:
% The transformed STM maps local state variations between orthogonal
% symplectic bases aligned with:
%
%   (i)   the local dynamical flow
%   (ii)  the local Hamiltonian gradient
%
% at the initial and final states of the trajectory segment.
%
% CURRENT IMPLEMENTATION:
%   - Planar 2-DoF subsystem only
%   - The planar subsystem is extracted internally from the full STM
%   - Full state vectors are passed directly to compute_R_mtx_2dof
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
% The full STM must correspond to the same formulation as the state vector:
%
%   Phi = d x_f / d x_0
%
% The planar subsystem:
%
%       idx = [1 2 4 5]
%
% is extracted internally.
%
% FORMULATION HANDLING:
% If the STM is provided in Lagrangian coordinates, it is internally
% transformed into canonical Hamiltonian coordinates before constructing
% Psi.
%
% INPUTS:                   Description                          Units
%
%   t0              - Initial time                              [-]
%   tf              - Final time                                [-]
%   x0_vec          - Initial full state vector (6x1)           [-]
%   xf_vec          - Final full state vector (6x1)             [-]
%   Phi_mtx         - Full STM (6x6)                            [-]
%   params          - Parameter struct                          [-]
%
% OUTPUTS:                  Description                          Units
%
%   Psi_mtx         - Transformed STM                           [-]
%
%
%==========================================================================

    % ---------------------------------------------------------------------
    %              EXTRACT PLANAR STM
    % ---------------------------------------------------------------------

    idx = [1, 2, 4, 5];

    Phi_planar_mtx = Phi_mtx(idx, idx);

    % ---------------------------------------------------------------------
    %              TRANSFORMATION TO CANONICAL COORDINATES
    % ---------------------------------------------------------------------

    % Transformation matrix: [x y vx vy] -> [x y px py]
    switch lower(params.model.name)
        case 'cr3bp'
            T_mtx = [eye(2), zeros(2);
                    [0,-1; 1, 0], eye(2)];
        case 'hillr3bp'
            T_mtx = [eye(2), zeros(2);
                    [0,-1; 1, 0], eye(2)];            
        case '2bp'
            T_mtx = eye(4);
        otherwise
            error('Unsupported model for compute_R_2dof.');
    end

    switch lower(params.model.formulation)

        case 'hamiltonian'

            Phi_ham_mtx = Phi_planar_mtx;

        case 'lagrangian'

            Phi_ham_mtx = T_mtx * Phi_planar_mtx / T_mtx;

        otherwise

            error('Unknown formulation');

    end

    % ---------------------------------------------------------------------
    %              COMPUTE ORTHO-SYMPLECTIC MATRICES
    % ---------------------------------------------------------------------

    R_o_mtx = compute_R_2dof(t0, x0_vec, params);

    R_f_mtx = compute_R_2dof(tf, xf_vec, params);

    % ---------------------------------------------------------------------
    %              COMPUTE TRANSFORMED STM
    % ---------------------------------------------------------------------

    Psi_mtx = R_f_mtx' * Phi_ham_mtx * R_o_mtx;

end