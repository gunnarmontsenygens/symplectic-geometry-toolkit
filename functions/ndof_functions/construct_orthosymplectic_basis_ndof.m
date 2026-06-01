function [u_hat_set, v_hat_set] = construct_orthosymplectic_basis_ndof(u1_tilde_vec)

    % Initialization
    u1_tilde_vec = u1_tilde_vec(:)';
    N = length(u1_tilde_vec);
    n = N/2;
    u_hat_set = zeros(n, N);
    v_hat_set = zeros(n, N);

    % Create the symplectic matrix
    J_mtx = [zeros(n), eye(n);
            -eye(n), zeros(n)];

    % Generate n linearly independent vectors
    u_tilde_vec_set = generate_LI_vec(u1_tilde_vec);

    % Assign the first u_hat
    u1_vec = u1_tilde_vec; u1_hat = u1_vec / norm(u1_vec);
    u_hat_set(1,:) = u1_hat;

    % Assign the j u_hat
    for j = 2 : n

        uj_tilde_vec = u_tilde_vec_set(j, :);

        % Calculate the terms in the sum first
        sum_j_vec = zeros(1,N);

        for k = 1 : j-1
            uk_hat = u_hat_set(k, :);
            sum_j_vec = sum_j_vec + dot(uj_tilde_vec, uk_hat)*uk_hat/norm(uk_hat)^2 ...
                + ((dot(uj_tilde_vec, J_mtx*uk_hat')/norm(J_mtx*uk_hat')^2)*J_mtx*uk_hat')';

        end
        
        % Append the new vectors
        uj_vec = uj_tilde_vec - sum_j_vec;
        uj_hat = uj_vec / norm(uj_vec);
        u_hat_set(j, :) = uj_hat;

    end

    % Compute v_hat_set
    for i = 1 : n

        v_hat_set(i, :) = - (J_mtx * u_hat_set(i, :)')';

    end
    
end