function [Phi_hist, PhiTPhi_hist, Psi_tilde_hist, Psi_tildeTPsi_tilde_hist, check_hist] = analyze_Phi_Psi_tilde_2dof(x0_vec, t_hist, params)
%==========================================================================
%
% Computes the time evolution of the State Transition Matrix (STM),
% transformed STM, stretching matrices, and associated symplectic structure
% diagnostics for a planar 2-DoF Hamiltonian subsystem.
%
% The transformed STM is defined as:
%
%   \tilde{\Psi}(t_f,t_0) = R_f^T \Phi(t_f,t_0) R_0
%
% where:
%
%   \Phi(t_f,t_0)    - State Transition Matrix (STM)
%   R_0              - Initial ortho-symplectic basis
%   R_f              - Final ortho-symplectic basis
%
% In addition to the STM and transformed STM, the function also computes:
%
%   \Phi^T \Phi
%   \tilde{\Psi}^T \tilde{\Psi}
%
% which characterize the Euclidean stretching induced by each map.
%
% MODEL DESCRIPTION:
% The transformed STM decouples the local phase space stretching into
% dynamically meaningful directions aligned with:
%
%   (i)   the local dynamical flow
%   (ii)  the local Hamiltonian gradient
%
% The function propagates the state and STM over the provided time history,
% extracts the planar STM, computes the transformed STM at each time step,
% evaluates the corresponding stretching matrices, and stores their
% eigenvalues and eigenvectors. It also evaluates structural and symplectic
% consistency checks associated with the transformed representation.
%
% CURRENT IMPLEMENTATION:
%   - Planar 2-DoF Hamiltonian subsystem only
%   - Full state vectors are propagated internally
%   - Planar STM extracted using:
%
%         idx = [1,2,4,5]
%
% OUTPUT STRUCTURES:
%
%   Phi_hist:
%       .mtx            - Planar STM history                    [-]
%       .evals          - STM eigenvalue history                [-]
%       .norm_evals     - STM eigenvalue norm history           [-]
%       .evecs          - STM eigenvector history               [-]
%
%   PhiTPhi_hist:
%       .mtx            - Phi^T Phi history                     [-]
%       .evals          - Phi^T Phi eigenvalue history          [-]
%       .norm_evals     - Phi^T Phi eigenvalue norms            [-]
%       .evecs          - Phi^T Phi eigenvector history         [-]
%
%   Psi_tilde_hist:
%       .mtx            - Transformed STM history               [-]
%       .evals          - Transformed STM eigenvalue history    [-]
%       .norm_evals     - Transformed STM eigenvalue norms      [-]
%       .evecs          - Transformed STM eigenvector history   [-]
%
%   Psi_tildeTPsi_tilde_hist:
%       .mtx            - Psi_tilde^T Psi_tilde history         [-]
%       .evals          - Psi_tilde^T Psi_tilde eigenvalues     [-]
%       .norm_evals     - Psi_tilde^T Psi_tilde eigenvalue norms[-]
%       .evecs          - Psi_tilde^T Psi_tilde eigenvectors    [-]
%
%   check_hist:
%       .sigma_vec              - Sigma vector history          [-]
%       .gamma_vec              - Gamma vector history          [-]
%       .lambda_vec             - Lambda vector history         [-]
%       .lambda_inv_check       - lambda*lambda^-1 check        [-]
%       .detM                   - det(M) history                [-]
%       .sympl_err              - Symplecticity error history   [-]
%       .zeroterms_err          - Zero-term constraint error    [-]
%       .sigmagamma_err_vec     - Sigma-gamma vector error      [-]
%       .sigmagamma_err         - Sigma-gamma norm error        [-]
%
% Author: G. Montseny
% Date: May 26, 2026
%
% INPUT:                    Description                          Units
%
%   x0_vec         - Initial state vector                       [-]
%   t_hist         - Time history                               [-]
%   params         - Parameter struct                           [-]
%
% OUTPUT:                   Description                          Units
%
%   Phi_hist       - STM diagnostic struct                      [-]
%   Psi_tilde_hist - Transformed STM diagnostic struct          [-]
%   check_hist     - Symplectic consistency diagnostic struct   [-]
%
%==========================================================================
    
    
    % Initialization
    x0_vec = x0_vec(:);
    n = length(x0_vec);
    N_t = length(t_hist);
    integrate = params.fun.integrate;
    
    % State and STM integration
    [~, x_vec_hist, Phi_mtx_hist, ~, ~, ~, ~] = integrate(t_hist, x0_vec, params);

    % Phi Preallocation
    Phi_hist = struct(...
        'mtx', zeros(N_t, 4, 4), ...
        'evals', zeros(N_t, 4), ...
        'norm_evals', zeros(N_t, 4), ...
        'evecs', zeros(N_t, 4, 4));

    %PhiTPhi Preallocation
    PhiTPhi_hist = struct(...
        'mtx', zeros(N_t, 4, 4), ...
        'evals', zeros(N_t, 4), ...
        'norm_evals', zeros(N_t, 4), ...
        'evecs', zeros(N_t, 4, 4));

    % Psi tilde Preallocation
    Psi_tilde_hist = struct(...
        'mtx', zeros(N_t, 4, 4), ...
        'evals', zeros(N_t, 4), ...
        'norm_evals', zeros(N_t, 4), ...
        'evecs', zeros(N_t, 4, 4));

    %PsiTildeTPsiTilde Preallocation
    Psi_tildeTPsi_tilde_hist = struct(...
        'mtx', zeros(N_t, 4, 4), ...
        'evals', zeros(N_t, 4), ...
        'norm_evals', zeros(N_t, 4), ...
        'evecs', zeros(N_t, 4, 4));

    % Check Psi Tilde Preallocation
    check_hist = struct( ...
        'sigma_vec', zeros(N_t,3), ...
        'gamma_vec', zeros(N_t,2), ...
        'lambda_vec', zeros(N_t,2), ...
        'lambda_inv_check', zeros(N_t,1), ...
        'detM', zeros(N_t,1), ...
        'sympl_err', zeros(N_t,1), ...
        'zeroterms_err', zeros(N_t,1), ...
        'sigmagamma_err_vec', zeros(N_t,2),...
        'sigmagamma_err', zeros(N_t,1));


    % TIME LOOP 
    for i = 1 : N_t

        % PHI -------------------------------------------------------------
        
        % Extract planar Phi
        Phi_mtx = squeeze(Phi_mtx_hist(i, :, :));
        idx = [1, 2, 4, 5];
        Phi_planar_mtx = Phi_mtx(idx, idx);

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

        % Calculate eigenvalues and eigenvectors of planar Phi
        [V,D] = eig(Phi_ham_mtx);
        
        % Store info in the respective arrays
        Phi_hist.mtx(i, :, :) = Phi_ham_mtx;
        Phi_hist.evals(i, :) = [D(1,1), D(2,2), D(3,3), D(4,4)];
        Phi_hist.norm_evals(i, :) = [norm(D(1,1)), norm(D(2,2)), norm(D(3,3)), norm(D(4,4))];
        Phi_hist.evecs(i, :, :) = V;

        % PHI^T PHI -------------------------------------------------------------
        
        % Compute Phi^T Phi
        PhiTPhi_mtx = Phi_ham_mtx' * Phi_ham_mtx;

        % Calculate eigenvalues and eigenvectors of planar Phi ^T Phi
        [V,D] = eig(PhiTPhi_mtx);
        
        % Store info in the respective arrays
        PhiTPhi_hist.mtx(i, :, :) = PhiTPhi_mtx;
        PhiTPhi_hist.evals(i, :) = [D(1,1), D(2,2), D(3,3), D(4,4)];
        PhiTPhi_hist.norm_evals(i, :) = [norm(D(1,1)), norm(D(2,2)), norm(D(3,3)), norm(D(4,4))];
        PhiTPhi_hist.evecs(i, :, :) = V;

        % PSI TILDE -------------------------------------------------------

        % Compute Psi_tilde matrix at t_hist(i)
        t0 = t_hist(1);
        tf = t_hist(i);
        xf_vec = x_vec_hist(i, :);
        Phi_f_mtx = squeeze(Phi_mtx_hist(i, :, :));
        Psi_tilde_mtx = compute_Psi_tilde_2dof(t0, tf, x0_vec, xf_vec, Phi_f_mtx, params);

        % Calculate eigenvalues and eigenvectors of Psi_tilde
        [V,D] = eig(Psi_tilde_mtx);        

        % Store info in the respective arrays
        Psi_tilde_hist.mtx(i, :, :) = Psi_tilde_mtx;
        Psi_tilde_hist.evals(i, :) = [D(1,1), D(2,2), D(3,3), D(4,4)];
        Psi_tilde_hist.norm_evals(i, :) = [norm(D(1,1)), norm(D(2,2)), norm(D(3,3)), norm(D(4,4))];
        Psi_tilde_hist.evecs(i, :, :) = V;

        % PSI TILDE ^T PSI TILDE -------------------------------------------------------

        % Compute 
        Psi_tildeTPsi_tilde_mtx = Psi_tilde_mtx' * Psi_tilde_mtx;

        % Calculate eigenvalues and eigenvectors
        [V,D] = eig(Psi_tildeTPsi_tilde_mtx);        

        % Store info in the respective arrays
        Psi_tildeTPsi_tilde_hist.mtx(i, :, :) = Psi_tildeTPsi_tilde_mtx;
        Psi_tildeTPsi_tilde_hist.evals(i, :) = [D(1,1), D(2,2), D(3,3), D(4,4)];
        Psi_tildeTPsi_tilde_hist.norm_evals(i, :) = [norm(D(1,1)), norm(D(2,2)), norm(D(3,3)), norm(D(4,4))];
        Psi_tildeTPsi_tilde_hist.evecs(i, :, :) = V;

        % CHECK PSI TILDE -------------------------------------------------
    
        % Compute check at t_hist(i)
        check = check_Psi_tilde_2dof(Psi_tilde_mtx, false);

        % Store 
        check_hist.sigma_vec(i,:) = check.sigma_vec.';
        check_hist.gamma_vec(i,:) = check.gamma_vec.';
        check_hist.lambda_vec(i,:) = check.lambda_vec.';
        check_hist.lambda_inv_check(i) = check.lambda_inv_check;
        check_hist.detM(i) = check.detM;
        check_hist.sympl_err(i) = check.sympl_err;
        check_hist.zeroterms_err(i) = check.zeroterms_err;
        check_hist.sigmagamma_err_vec(i,:) = check.sigmagamma_err_vec.';
        check_hist.sigmagamma_err(i) = check.sigmagamma_err;


    end
    


end