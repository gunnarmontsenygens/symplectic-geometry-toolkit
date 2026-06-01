function check = check_Psi_tilde_3dof(Psi_tilde_mtx, display)
%==========================================================================
%
% Checks the structural and symplectic properties of the transformed State
% Transition Matrix (STM) in coordinate-momentum pair ordering.
%
% MODEL DESCRIPTION:
% For a spatial 3-DoF Hamiltonian system, the transformed STM is written in
% the paired coordinate-momentum ordering:
%
%   [q_1 p_1 q_2 p_2 q_3 p_3]
%
% In this ordering, the expected structure is:
%
%   Psi_tilde =
%
%       [ lambda     sigma_0     sigma_1     sigma_2     sigma_3     sigma_4 ]
%       [   0       1/lambda        0           0           0           0    ]
%       [   0        gamma_1        .           .           .           .    ]
%       [   0        gamma_2        .           .           .           .    ]
%       [   0        gamma_3        .           .           .           .    ]
%       [   0        gamma_4        .           .           .           .    ]
%
% where the lower-right block is the reduced 4x4 matrix:
%
%   M = Psi_tilde(3:6,3:6)
%
% The expected constraints are:
%
%   lambda * (1/lambda) = 1
%
%   det(M) = 1
%
%   Psi_tilde^T * J_p * Psi_tilde = J_p
%
% together with the sigma-gamma coupling constraint:
%
%   sigma + lambda * M^T * J_4 * gamma = 0
%
% where:
%
%   sigma = [sigma_1 sigma_2 sigma_3 sigma_4]^T
%
%   gamma = [gamma_1 gamma_2 gamma_3 gamma_4]^T
%
% and:
%
%   J_4 =
%
%       [  0   1   0   0 ]
%       [ -1   0   0   0 ]
%       [  0   0   0   1 ]
%       [  0   0  -1   0 ]
%
% The symplectic structure is checked using the permuted symplectic
% identity matrix J_p corresponding to the paired ordering.
%
% Author: G. Montseny
% Date: May 29, 2026
%
% INPUTS:                   Description                          Units
%
%   Psi_tilde_mtx   - Transformed STM in paired ordering        [-]
%   display         - Display diagnostic output flag            [-]
%
% OUTPUTS:                  Description                          Units
%
%   check           - Struct containing diagnostic quantities    [-]
%
%       sigma_0             - Coupling scalar                    [-]
%       sigma_vec           - Coupling vector                    [-]
%       gamma_vec           - Coupling vector                    [-]
%       lambda_vec          - Stretching vector                  [-]
%       lambda_inv_check    - lambda*(1/lambda) - 1             [-]
%       detM_check          - det(M) - 1                         [-]
%       sympl_err           - Frobenius norm of symplectic defect[-]
%       zeroterms_err       - Sum of absolute zero-structure terms[-]
%       sigmagamma_err_vec  - Sigma-gamma coupling error vector  [-]
%       sigmagamma_err      - Norm of sigma-gamma coupling error [-]
%
%==========================================================================
    
    % If only one input is provided
    if nargin < 2
        display = true;
    end
    

    % Extract the parameters
    sigma_vec = [Psi_tilde_mtx(1,2:6)]';
    gamma_vec = [Psi_tilde_mtx(3:6,2)];
    lambda_vec = [Psi_tilde_mtx(1,1); Psi_tilde_mtx(2,2)];

    % Extract submatrix M
    M_mtx = Psi_tilde_mtx(3:6, 3:6);

    % Lambda * lambda ^-1
    lambda_inv_check = lambda_vec(1) * lambda_vec(2)-1;

    % detM
    detMm1 = det(M_mtx)-1;

    % Symplectic error using the frobenius norm
    J_p_mtx = compute_permuted_symplectic_identity_mtx_3dof();
    sympl_err = norm(Psi_tilde_mtx' * J_p_mtx * Psi_tilde_mtx - J_p_mtx, 'fro');


    % Sum of absolute values of entries that should be zero
    zeroterms_err = abs(Psi_tilde_mtx(2,1)) + abs(Psi_tilde_mtx(3,1)) + ...
                    abs(Psi_tilde_mtx(4,1)) + abs(Psi_tilde_mtx(5,1)) + ...
                    abs(Psi_tilde_mtx(6,1)) + abs(Psi_tilde_mtx(2,3)) + ...
                    abs(Psi_tilde_mtx(2,4)) + abs(Psi_tilde_mtx(2,5)) + ...
                    abs(Psi_tilde_mtx(2,6));

    % Checking the relationship between sigma_vec and gamma_vec
    J_4_mtx = [0,	1,	0,	0;
                -1,	0,	0,	0;
                0,	0,	0,	1;
                0,	0,	-1,	0];
    sigmagamma_err_vec = sigma_vec(2:5) + lambda_vec(1) * M_mtx.' * J_4_mtx * gamma_vec;
    sigmagamma_err = norm(sigmagamma_err_vec);

    % Incorporating all the values into the struct 
    check.sigma_vec = sigma_vec;
    check.gamma_vec = gamma_vec;
    check.lambda_vec = lambda_vec;
    check.lambda_inv_check = lambda_inv_check;
    check.detMm1 = detMm1;
    check.sympl_err = sympl_err;
    check.zeroterms_err = zeroterms_err;
    check.sigmagamma_err_vec = sigmagamma_err_vec;
    check.sigmagamma_err = sigmagamma_err;

    % Display checks
    if display
        disp(' ')
        disp('================ Psi Tilde Check ================')
        
        
        fprintf('sigma_vec = \n')
        fprintf('    %.15g\n', check.sigma_vec)
        
        fprintf('\ngamma_vec = \n')
        fprintf('    %.15g\n', check.gamma_vec)
        
        fprintf('\nlambda_vec = \n')
        fprintf('    %.15g\n', check.lambda_vec)
        
        fprintf('\nlambda_inv_check = %.15g\n', check.lambda_inv_check)
        fprintf('detMm1             = %.15g\n', check.detMm1)
        fprintf('sympl_err        = %.15g\n', check.sympl_err)
        fprintf('zeroterms_err    = %.15g\n', check.zeroterms_err)
        
        fprintf('\nsigmagamma_err_vec = \n')
        fprintf('    %.15g\n', check.sigmagamma_err_vec)
        
        fprintf('\nsigmagamma_err   = %.15g\n', check.sigmagamma_err)
        
        disp('=================================================')
    end

end