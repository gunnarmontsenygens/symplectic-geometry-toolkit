function Psi_mtx = compute_Psi_3dof(t0, tf, x0_vec, xf_vec, Phi_mtx, params)
%==========================================================================
%
% Computes the transformed State Transition Matrix (STM):
%
%   Psi = R_f^T * Phi * R_o
%
% for a spatial 3-DoF Hamiltonian dynamical system using the
% ortho-symplectic transformation introduced by Scheeres & Boodram (2025).
%
% MODEL DESCRIPTION:
% The transformed STM maps local state variations between orthogonal
% symplectic bases aligned with:
%
%   (i)    the local dynamical flow
%   (ii)   the local Hamiltonian gradient
%   (iii)  the remaining orthogonal symplectic directions
%
% at the initial and final states of the trajectory segment.
%
% The transformation is performed using the ortho-symplectic matrices:
%
%   R_o = R(t_0, x_0)
%   R_f = R(t_f, x_f)
%
% yielding:
%
%   Psi(t_f,t_0) = R_f^T Phi(t_f,t_0) R_o
%
% where:
%
%   Phi(t_f,t_0)   - State Transition Matrix (STM)
%   R_o            - Initial ortho-symplectic basis
%   R_f            - Final ortho-symplectic basis
%
% The resulting matrix expresses local phase-space stretching and rotation
% in a coordinate system naturally adapted to the Hamiltonian flow.
%
% CURRENT IMPLEMENTATION:
%   - Spatial 3-DoF systems only
%   - Full 6-dimensional state vectors are used directly
%   - Supports both Lagrangian and Hamiltonian formulations
%   - Lagrangian STMs are internally transformed to canonical coordinates
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
% where:
%
%       x_0 = x(t_0)
%       x_f = x(t_f)
%
% FORMULATION HANDLING:
% If the STM is provided in Lagrangian coordinates, it is internally
% transformed into canonical Hamiltonian coordinates prior to computing
% Psi:
%
%       Phi_ham = T Phi_lag T^{-1}
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
%   Psi_mtx         - Transformed STM                           [-]
%
%==========================================================================


    % ---------------------------------------------------------------------
    %              TRANSFORMATION TO CANONICAL COORDINATES
    % ---------------------------------------------------------------------

    % Transformation matrix: [x y vx vy] -> [x y px py]
    switch lower(params.model.name)
        case 'cr3bp'
            T_mtx = [eye(3), zeros(3);
                    [0,-1, 0; 1, 0, 0; 0, 0, 0], eye(3)];
        case 'hillr3bp'
            T_mtx = [eye(3), zeros(3);
                    [0,-1, 0; 1, 0, 0; 0, 0, 0], eye(3)];  
        case 'apccr4bp'
            T_mtx = [eye(3), zeros(3);
                    [0,-1, 0; 1, 0, 0; 0, 0, 0], eye(3)]; 
        case '2bp'
            T_mtx = eye(6);
        otherwise
            error('Unsupported model for compute_R_3dof.');
    end

    switch lower(params.model.formulation)

        case 'hamiltonian'

            Phi_ham_mtx = Phi_mtx;

        case 'lagrangian'

            Phi_ham_mtx = T_mtx * Phi_mtx / T_mtx;

        otherwise

            error('Unknown formulation');

    end

    % ---------------------------------------------------------------------
    %              COMPUTE ORTHO-SYMPLECTIC MATRICES
    % ---------------------------------------------------------------------

    R_o_mtx = compute_R_3dof(t0, x0_vec, params);

    R_f_mtx = compute_R_3dof(tf, xf_vec, params);

    % ---------------------------------------------------------------------
    %              COMPUTE TRANSFORMED STM
    % ---------------------------------------------------------------------

    Psi_mtx = R_f_mtx' * Phi_ham_mtx * R_o_mtx;

end