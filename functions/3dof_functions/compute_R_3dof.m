function R_mtx = compute_R_3dof(t, x_vec, params)
%==========================================================================
%
% Computes the ortho-symplectic transformation matrix for a spatial
% 3-DoF Hamiltonian dynamical system using the closed-form construction
% introduced by Scheeres & Boodram (2025).
%
% MODEL DESCRIPTION:
% The transformation is constructed from the local Hamiltonian flow and
% Hamiltonian gradient evaluated at a given state. The resulting matrix
% separates the phase-space directions associated with:
%
%   (i)    the local dynamical flow
%   (ii)   the local Hamiltonian gradient
%   (iii)  the remaining orthogonal symplectic directions
%
% For a 3-DoF system, the transformation matrix is:
%
%   R = [u_1  u_2  u_3  v_1  v_2  v_3]
%
% where the basis vectors form an orthonormal symplectic basis satisfying:
%
%   R^T R = I
%   R^T J R = J
%
% The first vector pair is aligned with the Hamiltonian gradient and flow:
%
%   v_1 = H_x / ||H_x||
%   u_1 = J v_1
%
% The remaining vector pairs are constructed from the configuration-space
% and momentum-space components of the Hamiltonian gradient using the normal
% direction:
%
%   n_pi = H_q x H_p
%
% CURRENT IMPLEMENTATION:
%   - Spatial 3-DoF systems only
%   - Uses the explicit closed-form construction from the reference paper
%   - Supports full 6-dimensional state vectors
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
% FORMULATION HANDLING:
% If the system is provided in Lagrangian coordinates, the local flow is
% transformed internally into canonical Hamiltonian coordinates before
% constructing the ortho-symplectic basis.
%
% Author: G. Montseny
% Date: May 29, 2026
%
% INPUTS:                   Description                          Units
%
%   t               - Time                                      [-]
%   x_vec           - Full state vector (6x1)                   [-]
%   params          - Parameter struct                          [-]
%
% OUTPUTS:                  Description                          Units
%
%   R_mtx           - Ortho-symplectic transformation matrix     [-]
%
%==========================================================================

    % ---------------------------------------------------------------------
    %              COMPUTE FLOW
    % ---------------------------------------------------------------------

    % Use Model EoM
    dx_dt_vec = params.fun.eom(t, x_vec, params);

    % ---------------------------------------------------------------------
    %              TRANSFORMATION TO CANONICAL COORDINATES
    % ---------------------------------------------------------------------
    
    % Transformation matrix
    switch lower(params.model.name)
        case 'cr3bp'
            T_mtx = [eye(3), zeros(3);
                    [0,-1, 0; 1, 0, 0], eye(3)];
        case 'hillr3bp'
            T_mtx = [eye(3), zeros(3);
                    [0,-1, 0; 1, 0, 0], eye(3)];   
        case '2bp'
            T_mtx = eye(6);
        otherwise
            error('Unsupported model for compute_R_3dof.');
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
        
    J = [zeros(3,3), eye(3);
        - eye(3), zeros(3,3)];
    H_x_vec = - J * dx_dt_vec;

    % ---------------------------------------------------------------------
    % (ii) Calculate unit vector and essential components
    % ---------------------------------------------------------------------
    
    H_q_vec = H_x_vec(1:3);
    H_p_vec = H_x_vec(4:6);  
    H_x = norm(H_x_vec);
    H_x_hat = H_x_vec/H_x; 
    
    % ---------------------------------------------------------------------
    % (iii) Calculate u_1 and v_1
    % ---------------------------------------------------------------------
    
    u_1_hat = J * H_x_hat;
    v_1_hat = H_x_hat;

    % ---------------------------------------------------------------------
    % (iv) Calculate u_3 and v_3
    % ---------------------------------------------------------------------
    
    % Calculate normal to the plane
    n_pi_vec = cross(H_q_vec, H_p_vec);
    n_pi_hat = n_pi_vec / norm(n_pi_vec);

    % Calculate u_3 and v_3
    u_3_hat = [n_pi_hat; zeros(3,1)];
    v_3_hat = [zeros(3,1); n_pi_hat];

    % ---------------------------------------------------------------------
    % (iv) Calculate u_2 and v_2
    % ---------------------------------------------------------------------
    
    u_2_hat = [cross(n_pi_hat, H_p_vec); cross(n_pi_hat, H_q_vec)] / H_x;
    v_2_hat = [ - cross(n_pi_hat, H_q_vec); cross(n_pi_hat, H_p_vec)] / H_x;

    % ---------------------------------------------------------------------
    % (v) Calculate R_mtx
    % ---------------------------------------------------------------------
    R_mtx = [u_1_hat, u_2_hat, u_3_hat v_1_hat, v_2_hat, v_3_hat];

end