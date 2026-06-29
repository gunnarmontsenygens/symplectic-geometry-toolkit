function [Psi_tilde_mtx, Psi_tilde_red_mtx] = compute_Psi_tilde_moreno_2dof(t0, tf, x0_vec, xf_vec, Phi_mtx, params)

    % ---------------------------------------------------------------------
    %              PRELIMINARIES
    % ---------------------------------------------------------------------

    % Indexing
    switch lower(params.model.name)
        case 'pccr4bp'
            idx = [1,2,3,4];
        otherwise
            idx = [1, 2, 4, 5];
    end

    % Extract planar STM
    Phi_planar_mtx = Phi_mtx(idx, idx);
    
    % Transform into canonical coordinates: [x y vx vy] -> [x y px py]
    switch lower(params.model.name)
        case 'cr3bp'
            T_mtx = [eye(2), zeros(2);
                    [0,-1; 1, 0], eye(2)];
        case 'hillr3bp'
            T_mtx = [eye(2), zeros(2);
                    [0,-1; 1, 0], eye(2)];        
        case 'pccr4bp'
            T_mtx = [eye(2), zeros(2);
                    [0,-1; 1, 0], eye(2)];       
        case '2bp'
            T_mtx = eye(4);
        otherwise
            error('Unsupported model for compute_R_2dof.');
    end

    % Compute hamiltonian version of the STM
    switch lower(params.model.formulation)
        case 'hamiltonian'
            Phi_ham_mtx = Phi_planar_mtx;

        case 'lagrangian'
            Phi_ham_mtx = T_mtx * Phi_planar_mtx / T_mtx;

        otherwise
            error('Unknown formulation');
    end


    % ---------------------------------------------------------------------
    %              COMPUTE F & F BAR
    % ---------------------------------------------------------------------

    % Compute F and F_bar at t = t0
    [F0_mtx, F0_bar_mtx] = compute_FF_bar_2dof(t0, x0_vec, params)


    % Compute F and F_bar at t = tf
    [Ff_mtx, Ff_bar_mtx] = compute_FF_bar_2dof(tf, xf_vec, params)


    % ---------------------------------------------------------------------
    %              COMPUTE PSI and PSI_RED
    % ---------------------------------------------------------------------

    % Define identity symplectic matrices
    Omega_4_mtx = [zeros(2,2), eye(2);
            - eye(2), zeros(2,2)];
    Omega_2_mtx = [zeros(1,1), eye(1);
            - eye(1), zeros(1,1)];

    % Compute Psi_mtx
    Psi_mtx = Omega_4_mtx.' * Ff_bar_mtx.' * Omega_4_mtx * Phi_ham_mtx * F0_bar_mtx;

    % Compute Psi_red_mtx 
    Psi_red_mtx = Omega_2_mtx.' * Ff_mtx.' * Omega_4_mtx * Phi_ham_mtx * F0_mtx;

    % ---------------------------------------------------------------------
    %              COMPUTE P_x2z
    % ---------------------------------------------------------------------

    P_x2z_mtx = compute_P_x2z_2dof();

    % ---------------------------------------------------------------------
    %              COMPUTE PSI_TILDE & PSI_TILDE_RED
    % ---------------------------------------------------------------------

    Psi_tilde_mtx = P_x2z_mtx * Psi_mtx * P_x2z_mtx';

    Psi_tilde_red_mtx =   Psi_red_mtx ;

    % ---------------------------------------------------------------------
    %              HELPER FUNCTIONS
    % ---------------------------------------------------------------------

    function [F_mtx, F_bar_mtx] = compute_FF_bar_2dof(t, x_vec, params)

        % Indexing
            switch lower(params.model.name)
                case 'pccr4bp'
                    idx_F = [1,2,3,4];
                otherwise
                    idx_F = [1, 2, 4, 5];
            end 

        % Transform into canonical coordinates: [x y vx vy] -> [x y px py]
        switch lower(params.model.name)
            case 'cr3bp'
                T_F_mtx = [eye(2), zeros(2);
                        [0,-1; 1, 0], eye(2)];
            case 'hillr3bp'
                T_F_mtx = [eye(2), zeros(2);
                        [0,-1; 1, 0], eye(2)];        
            case 'pccr4bp'
                T_F_mtx = [eye(2), zeros(2);
                        [0,-1; 1, 0], eye(2)];       
            case '2bp'
                T_F_mtx = eye(4);
            otherwise
                error('Unsupported model for compute_R_2dof.');
        end
        
        % Extract dx_dt
        dx_dt_vec = params.fun.eom(t, x_vec, params);
        dx_dt_vec = dx_dt_vec(idx_F);
    
        switch params.model.formulation
            case 'hamiltonian'
            case 'lagrangian'
                dx_dt_vec = T_F_mtx * dx_dt_vec;
    
            otherwise
                error('Unknown formulation');
        end

        % (i) Calculate the first two components
        X_H_vec = dx_dt_vec; X_H_vec = X_H_vec(:);
    
        Omega_4 = [zeros(2,2), eye(2);
            - eye(2), zeros(2,2)];
        gradH_vec = -Omega_4 * dx_dt_vec;
    
        Z_vec = gradH_vec / norm(gradH_vec)^2; Z_vec = Z_vec(:);
    
        % (ii) Create quaternionic matrices
        I_mtx = [0, 0, -1, 0;
                0, 0, 0, -1;
                1, 0, 0, 0;
                0, 1, 0, 0];
    
        J_mtx = [0, -1, 0, 0;
                1, 0, 0, 0;
                0, 0, 0, 1;
                0, 0, -1, 0];
    
        K_mtx = [0, 0, 0, -1;
                0, 0, 1, 0;
                0, -1, 0, 0;
                1, 0, 0, 0];
    
        % (iii) Create the last two components
        U1_vec =  - J_mtx * gradH_vec / norm(gradH_vec)^2; U1_vec = U1_vec(:);
        V1_vec = K_mtx * gradH_vec; V1_vec = V1_vec(:);
    
        % (iv) Compute F
        F_mtx = [U1_vec, V1_vec];
    
        % (v) Compute F_bar
        F_bar_mtx = [Z_vec, U1_vec, X_H_vec, V1_vec];

    end


end