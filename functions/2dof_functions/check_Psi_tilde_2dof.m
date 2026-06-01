function check = check_Psi_tilde_2dof(Psi_tilde_mtx, display)
%==========================================================================
%
% Checks the structural and symplectic properties of the transformed State
% Transition Matrix (STM) in the coordinate-momentum pair ordering.
%
% MODEL DESCRIPTION:
% For a planar 2-DoF Hamiltonian subsystem, the transformed STM has the
% expected structure:
%
%   Psi_tilde =
%
%       [ lambda     sigma_0   sigma_1   sigma_2 ]
%       [   0       1/lambda      0         0    ]
%       [   0        gamma_1      a         b    ]
%       [   0        gamma_2      c         d    ]
%
% where:
%
%   M = [a b; c d]
%
% is the reduced 2x2 symplectic submatrix.
%
% The expected constraints are:
%
%   lambda * (1/lambda) = 1
%
%   det(M) = 1
%
%   Psi_tilde^T * J * Psi_tilde = J
%
%   gamma = M * J_2 * sigma
%
% where:
%
%   sigma = [sigma_1; sigma_2]
%
%   gamma = [gamma_1; gamma_2]
%
%   J_2 =
%
%       [  0   1 ]
%       [ -1   0 ]
%
% Author: G. Montseny
% Date: May 14, 2026
%
% INPUTS:                   Description                          Units
%
%   Psi_tilde_mtx  - Transformed STM in paired ordering         [-]
%
% OUTPUTS:                  Description                          Units
%
%   sigma_vec      - Coupling vector [sigma_0; sigma_1; sigma_2] [-]
%   gamma_vec      - Coupling vector [gamma_1; gamma_2]          [-]
%   lambda_vec     - Stretching vector [lambda; 1/lambda]        [-]
%   lambda_inv_check - Product lambda * (1/lambda)               [-]
%   detM           - Determinant of reduced matrix M             [-]
%   sympl_err      - Frobenius norm of symplectic defect         [-]
%   zeroterms_err  - Sum of absolute zero-structure entries      [-]
%   sigmagamma_err_vec - Error in gamma = M * J_2 * sigma        [-]
%   sigmagamma_err - Norm of error in gamma = M * J_2 * sigma    [-]
%
%
%==========================================================================
    
    % If only one input is provided
    if nargin < 2
        display = true;
    end
    

    % Extract the parameters
    sigma_vec = [Psi_tilde_mtx(1,2); Psi_tilde_mtx(1,3); Psi_tilde_mtx(1,4)];
    gamma_vec = [Psi_tilde_mtx(3,2); Psi_tilde_mtx(4,2)];
    lambda_vec = [Psi_tilde_mtx(1,1); Psi_tilde_mtx(2,2)];

    % Extract submatrix M
    M_mtx = Psi_tilde_mtx(3:4, 3:4);

    % Lambda * lambda ^-1
    lambda_inv_check = lambda_vec(1) * lambda_vec(2)-1;

    % detM
    detM = det(M_mtx)-1;

    % Symplectic error using the frobenius norm
    J_mtx = [ 0,  1,  0,  0;
       -1,  0,  0,  0;
        0,  0,  0,  1;
        0,  0, -1,  0];

    sympl_err = norm(Psi_tilde_mtx' * J_mtx * Psi_tilde_mtx - J_mtx, 'fro');


    % Sum of absolute values of entries that should be zero
    zeroterms_err = abs(Psi_tilde_mtx(2,1)) + abs(Psi_tilde_mtx(3,1)) + ...
                    abs(Psi_tilde_mtx(4,1)) + abs(Psi_tilde_mtx(2,3)) + ...
                    abs(Psi_tilde_mtx(2,4));

    % Checking the relationship between sigma_vec and gamma_vec
    J_2_mtx = [0, 1;
            -1, 0];
    sigmagamma_err_vec = gamma_vec(1:2) - Psi_tilde_mtx(2,2)*M_mtx * J_2_mtx * sigma_vec(2:3);
    sigmagamma_err = norm(sigmagamma_err_vec);

    % Incorporating all the values into the struct 
    check.sigma_vec = sigma_vec;
    check.gamma_vec = gamma_vec;
    check.lambda_vec = lambda_vec;
    check.lambda_inv_check = lambda_inv_check;
    check.detM = detM;
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
        fprintf('detM             = %.15g\n', check.detM)
        fprintf('sympl_err        = %.15g\n', check.sympl_err)
        fprintf('zeroterms_err    = %.15g\n', check.zeroterms_err)
        
        fprintf('\nsigmagamma_err_vec = \n')
        fprintf('    %.15g\n', check.sigmagamma_err_vec)
        
        fprintf('\nsigmagamma_err   = %.15g\n', check.sigmagamma_err)
        
        disp('=================================================')
    end

end