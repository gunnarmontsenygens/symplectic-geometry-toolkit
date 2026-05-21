function R_mtx = compute_R_ndof(t, x_vec, params)

    % ---------------------------------------------------------------------
    %              INITIALIZATION
    % ---------------------------------------------------------------------
    x_vec = x_vec(:);
    N = length(x_vec);
    n = N / 2;
    
    % ---------------------------------------------------------------------
    %              COMPUTE FLOW
    % ---------------------------------------------------------------------

    % Use Model EoM
    dx_dt_vec = params.fun.eom(t, x_vec, params);
    dx_dt_vec = dx_dt_vec(:);

    % ---------------------------------------------------------------------
    %              TRANSFORMATION TO CANONICAL COORDINATES
    % ---------------------------------------------------------------------
    
    % Transformation matrix
    switch lower(params.model.name)
        case 'cr3bp'
            if n == 2
                T_mtx = [eye(n), zeros(n);
                    [0,-1; 1, 0 ], eye(n)];                
            elseif n == 3
                 T_mtx = [eye(n), zeros(n);
                    [0,-1, 0; 1, 0, 0], eye(n)];                
            end
        case 'hillr3bp'
            if n == 2
                T_mtx = [eye(n), zeros(n);
                    [0,-1; 1, 0 ], eye(n)];                
            elseif n == 3
                 T_mtx = [eye(n), zeros(n);
                    [0,-1, 0; 1, 0, 0], eye(n)];                
            end
        case '2bp'
            T_mtx = eye(n);
        otherwise
            error('Unsupported model for compute_R_ndof.');
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
        
    J_mtx = compute_symplectic_identity_mtx(n);
    H_x_vec = - J_mtx * dx_dt_vec;
    H_x_hat = H_x_vec / norm(H_x_vec);

    % ---------------------------------------------------------------------
    % (ii) Calculate u1_tilde_vec
    % ---------------------------------------------------------------------
    
    u1_tilde_vec = J * H_x_hat;
    
    % ---------------------------------------------------------------------
    % (iii) Calculate u_hat_set and v_hat_set
    % ---------------------------------------------------------------------
    
    [u_hat_set, v_hat_set] = construct_orthosymplectic_basis_ndof(u1_tilde_vec);

    % ---------------------------------------------------------------------
    % (v) Calculate R_mtx
    % ---------------------------------------------------------------------
    
    R_mtx = zeros(N, N);

    for i = 1 : n
        R_mtx(:, i) = u_hat_set(i,:)';
        R_mtx(:, n + i) = v_hat_set(i,:)';
    end
    
end