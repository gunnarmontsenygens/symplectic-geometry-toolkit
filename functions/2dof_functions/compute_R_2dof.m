function R_mtx = compute_R_2dof(t, x_vec, params)
%==========================================================================
%
% Computes the ortho-symplectic transformation matrix for a planar
% 2-DoF Hamiltonian dynamical system using the closed-form construction
% introduced by Scheeres & Boodram (2025).
%
% MODEL DESCRIPTION:
% The transformation is constructed from the local Hamiltonian flow and
% Hamiltonian gradient evaluated at a given state. The resulting matrix
% separates the phase-space directions associated with:
%
%   (i)   the local dynamical flow
%   (ii)  the local Hamiltonian gradient
%
% from the remaining orthogonal symplectic directions.
%
% The transformation matrix:
%
%   R = [u_1  u_2  v_1  v_2]
%
% is both orthogonal and symplectic, satisfying:
%
%   R^T R = I
%   R^T J R = J
%
% CURRENT IMPLEMENTATION:
%   - Planar 2-DoF subsystem only
%   - Uses the explicit closed-form construction from the reference paper
%   - The planar subsystem is extracted internally from the full state
%
% STATE DEFINITION:
% The input state vector must be the full 6-dimensional state vector of the
% dynamical model.
%
% Supported state formulations:
%
%   Lagrangian:
%
%       x_vec = [x y z vx vy vz]^T
%
%   Hamiltonian:
%
%       x_vec = [x y z px py pz]^T
%
% The planar 2-DoF subsystem:
%
%       [x y vx vy]
%
% or
%
%       [x y px py]
%
% is extracted internally using:
%
%       idx = [1 2 4 5]
%
% FORMULATION HANDLING:
% If the system is provided in Lagrangian coordinates, the local flow is
% transformed internally into canonical Hamiltonian coordinates before
% constructing the ortho-symplectic basis.
%
% Author: G. Montseny
% Date: May 12, 2026
%
% INPUTS:                   Description                          Units
%
%   t               - Time                                      [-]
%   x_vec           - Full state vector (6x1)                   [-]
%   params          - Parameter struct                          [-]
%
% OUTPUTS:                  Description                          Units
%
%   R_mtx               - Ortho-symplectic transformation matrix    [-]
%
%==========================================================================

    % ---------------------------------------------------------------------
    %              COMPUTE FLOW
    % ---------------------------------------------------------------------

    % Use Model EoM
    dx_dt_vec = params.fun.eom(t, x_vec, params);

    % 2 DoF case
    switch lower(params.model.name)
        case 'pccr4bp'
            idx = [1,2,3,4];
        otherwise
            idx = [1, 2, 4, 5];
    end

    dx_dt_vec = dx_dt_vec(idx);

    % ---------------------------------------------------------------------
    %              TRANSFORMATION TO CANONICAL COORDINATES
    % ---------------------------------------------------------------------
    
    % Transformation matrix
    switch lower(params.model.name)
        case 'cr3bp'
            T_mtx = [eye(2), zeros(2);
                    [0,-1; 1, 0], eye(2)];

        case 'pccr4bp'
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


    % Transformation to canonical coordinates
    switch params.model.formulation
        case 'hamiltonian'
        case 'lagrangian'
            dx_dt_vec = T_mtx*dx_dt_vec;
        otherwise
            error('Unknown formulation');
    end

    % ---------------------------------------------------------------------
    % (i) Extract the Hamiltonian gradient
    % ---------------------------------------------------------------------
        
    J = [zeros(2,2), eye(2);
        - eye(2), zeros(2,2)];
    H_x_vec = -J*dx_dt_vec;

    % ---------------------------------------------------------------------
    % (ii) Calculate unit vector and essential components
    % ---------------------------------------------------------------------
    
    H_q_vec = H_x_vec(1:2);
    H_p_vec = H_x_vec(3:4);  
    H_x = norm(H_x_vec);
    H_x_hat = H_x_vec/H_x; 
    
    % ---------------------------------------------------------------------
    % (iii) Calculate u_1 and v_1
    % ---------------------------------------------------------------------
    
    u_1_hat = J*H_x_hat;
    v_1_hat = H_x_hat;

    % ---------------------------------------------------------------------
    % (iv) Calculate u_2 and v_2
    % ---------------------------------------------------------------------

    z_hat = [0; 0; 1];
    H_q_vec = [H_q_vec; 0];
    H_p_vec = [H_p_vec; 0];
    zcrossHp_vec = cross(z_hat,H_p_vec);
    zcrossHq_vec = cross(z_hat,H_q_vec);
    u_2_vec = [zcrossHp_vec(1:2); zcrossHq_vec(1:2)];
    u_2_hat = u_2_vec/H_x;
    v_2_hat = -J*u_2_hat;

    % ---------------------------------------------------------------------
    % (v) Calculate R_mtx
    % ---------------------------------------------------------------------
    R_mtx = [u_1_hat, u_2_hat, v_1_hat, v_2_hat];
end