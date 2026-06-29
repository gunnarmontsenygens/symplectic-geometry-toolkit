function P_x2z_mtx = compute_P_x2z_2dof()
%==========================================================================
%
% Computes the permutation matrix that transforms the canonical state
% ordering:
%
%   x = [q_1 q_2 p_1 p_2]^T
%
% into the coordinate-momentum pair ordering:
%
%   z = [q_1 p_1 q_2 p_2]^T
%
% MODEL DESCRIPTION:
% The ortho-symplectic transformation introduced by Scheeres & Boodram
% (2025) is naturally constructed in the canonical ordering:
%
%   [q_1 q_2 p_1 p_2]
%
% However, the transformed STM can be more conveniently analyzed in the
% pairwise coordinate-momentum ordering:
%
%   [q_1 p_1 q_2 p_2]
%
% This permutation matrix enables the transformation:
%
%   z = P_x2z * x
%
% and is used to construct:
%
%   Psi_tilde = P_x2z * Psi * P_x2z^T
%
% CURRENT IMPLEMENTATION:
%   - Planar 2-DoF systems only
%
% Author: G. Montseny
% Date: May 12, 2026
%
% INPUTS:                   Description                          Units
%
%   None
%
% OUTPUTS:                  Description                          Units
%
%   P_x2z_mtx          - Permutation matrix                       [-]
%
%
%==========================================================================
    
    P_x2z_mtx = [1, 0, 0, 0;
                0, 0, 1, 0;
                0, 1, 0, 0;
                0, 0, 0, 1];

end