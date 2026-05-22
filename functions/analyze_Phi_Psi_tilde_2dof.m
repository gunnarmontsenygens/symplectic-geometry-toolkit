function [Phi_mtx_hist, evals_Phi_mtx_hist, norm_evals_Phi_mtx_hist, evecs_Phi_mtx_hist, Psi_tilde_mtx_hist, evals_Psi_tilde_mtx_hist, norm_evals_Psi_tilde_mtx_hist, evecs_Psi_tilde_mtx_hist, check_hist] = analyze_Phi_Psi_tilde_2dof(x0_vec, t_hist, params)
    
    % ADD TOP COMMENT AND MAYBE GIVE STRUCT ARRAYS AS ANSWERS
    
    % Initialization
    x0_vec = x0_vec(:);
    n = length(x_0_vec);
    N_t = length(t_hist);
    integrate = params.fun.integrate;
    
    % State and STM integration
    [~, x_vec_hist, Phi_mtx_hist, ~, ~, ~, ~] = integrate(t_hist, x0_vec, params);

    % Phi Preallocation
    evals_Phi_mtx_hist = zeros(N_t, 4);
    norm_evals_Phi_mtx_hist = zeros(N_t, 4);
    evecs_Phi_mtx_hist = zeros(N_t, 4, 4);

    % Psi tilde Preallocation
    Psi_tilde_mtx_hist = zeros(N_t, 4, 4);
    evals_Psi_tilde_mtx_hist = zeros(N_t, 4);
    norm_evals_Psi_tilde_mtx_hist = zeros(N_t, 4);
    evecs_Psi_tilde_mtx_hist = zeros(N_t, 4, 4);

    % Check Psi Tilde Preallocation
    check_hist(N_t) = struct( ...
        'sigma_vec', [], ...
        'gamma_vec', [], ...
        'lambda_vec', [], ...
        'lambda_inv_check', [], ...
        'detM', [], ...
        'sympl_err', [], ...
        'zeroterms_err', [], ...
        'sigmagamma_err_vec', ...
        'sigmagamma_err', []);


    % TIME LOOP 
    for i = 1 : N_t

        % PHI -------------------------------------------------------------
        
        % Extract planar Phi
        Phi_mtx = squeeze(Phi_mtx_hist(N_t, :, :));
        idx = [1, 2, 4, 5];
        Phi_planar_mtx = Phi_mtx(idx, idx);

        % Calculate eigenvalues and eigenvectors of planar Phi
        [V,D] = eig(Phi_planar_mtx);
        
        % Store info in the respective arrays
        evals_Phi_mtx_hist(i, :) = [D(1,1), D(2,2), D(3,3), D(4,4)];
        norm_evals_Phi_mtx_hist(i, :) = [norm(D(1,1)), norm(D(2,2)), norm(D(3,3)), norm(D(4,4))];
        evecs_Phi_mtx_hist(i, :, :) = V;

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
        Psi_tilde_mtx_hist(i, :, :) = Psi_tilde_mtx;
        evals_Psi_tilde_mtx_hist(i, :) = [D(1,1), D(2,2), D(3,3), D(4,4)];
        norm_evals_Psi_tilde_mtx_hist(i, :) = [norm(D(1,1)), norm(D(2,2)), norm(D(3,3)), norm(D(4,4))];
        evecs_Psi_tilde_mtx_hist(i, :, :) = V;

        % CHECK PSI TILDE -------------------------------------------------
    
        % Compute check at t_hist(i)
        check = check_Psi_tilde_2dof(Psi_tilde_mtx);

        % Store 
        check_hist(i).sigma_vec = check.sigma_vec;
        check_hist(i).gamma_vec = check.gamma_vec;
        check_hist(i).lambda_vec = check.lambda_vec;
        check_hist(i).lambda_inv_check = check.lambda_inv_check;
        check_hist(i).detM = check.detM;
        check_hist(i).sympl_err = check.sympl_err;
        check_hist(i).zeroterms_err = check.zeroterms_err;
        check_hist(i).sigmagamma_err_vec = check.sigmagamma_err_vec;
        check_hist(i).sigmagamma_err = check.sigmagamma_err;


    end
    


end