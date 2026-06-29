function winding = compute_winding_number(S_mtx_hist, unitary)
%==========================================================================
%
% Computes phase-based winding diagnostics associated with a time history
% of symplectic matrices.
%
% Depending on the value of the input flag `unitary`, the computation is
% performed using either:
%
%   unitary = true
%
%       The orthogonal-symplectic factor obtained from the polar
%       decomposition
%
%           S = U P,
%
%       where U is extracted via the Singular Value Decomposition (SVD)
%
%           S = V_1 D V_2^T,
%           U = V_1 V_2^T.
%
%   unitary = false
%
%       The original symplectic matrix S itself.
%
% In either case, the matrix is first transformed to the standard
% (q,p)-ordered symplectic basis and written in block form
%
%       [ A  -B ]
%       [ B  A ].
%
% The associated complex matrix is then defined as
%
%       C = A + iB.
%
% For orthogonal-symplectic matrices this corresponds to the standard
% identification
%
%       U(n) ≅ Sp(2n,R) ∩ O(2n).
%
% The eigenvalues of the complex matrix are computed:
%
%       lambda_i = rho_i exp(i theta_i),
%
% together with their phase angles
%
%       theta_i = Im(log(lambda_i)).
%
% The determinant
%
%       d = det(C)
%
% provides a global phase whose argument equals the sum of the individual
% eigenvalue phases:
%
%       arg(d) = sum_i theta_i.
%
% The determinant phase history is unwrapped in time and the total winding
% angle is computed as
%
%       theta_tot_end = theta_tot_hist(t_f) - theta_tot_hist(t_0).
%
% unitary = true: winding of the orthogonal-symplectic (unitary) factor.
% unitary = false: evolution of the eigenvalue arguments of the full symplectic matrix S.
%
% OUTPUT STRUCTURE:
%
%   winding:
%
%       .theta_tot_hist
%           Unwrapped determinant phase history               [deg]
%
%       .theta_vec_hist
%           Unwrapped eigenvalue phase histories              [deg]
%
%       .lambda_vec_hist
%           Eigenvalue histories of C                         [-]
%
%       .C_mtx_hist
%           Complex matrix histories                          [-]
%
%       .theta_tot_end
%           Total winding angle                               [deg]
%
%       .unit_err_hist
%           ||C^* C - I||_F history                           [-]
%
%       .det_err_hist
%           |||det(C)| - 1| history                           [-]
%
% Author: G. Montseny
% Date: June 24, 2026
%
% INPUT:                    Description                          Units
%
%   S_mtx_hist     - Symplectic matrix history (N_t x 2n x 2n)   [-]
%   unitary        - Use orthogonal-symplectic factor from the
%                    polar decomposition if true; otherwise use
%                    the original symplectic matrix              [-]
%
% OUTPUT:                   Description                          Units
%
%   winding        - Winding diagnostic structure               [-]
%
%==========================================================================

    % Initialization
    N_t = size(S_mtx_hist,1);
    n = size(S_mtx_hist,2)/2;

    theta_hist = zeros(N_t,1);
    unit_err_hist = zeros(N_t,1);
    det_err_hist = zeros(N_t,1);
    if unitary
        lambda_vec_hist = zeros(N_t, n);
        theta_vec_hist = zeros(N_t, n);
        C_mtx_hist = zeros(N_t, n, n);
    else
        lambda_vec_hist = zeros(N_t, 2*n);
        theta_vec_hist = zeros(N_t, 2*n);
        C_mtx_hist = zeros(N_t, 2*n, 2*n);
    end

    I_n = eye(n);
    I_2n = eye(2*n);

    % Loop over history
    for k = 1:N_t

        % Extract symplectic matrix
        S_mtx = squeeze(S_mtx_hist(k,:,:));

        if unitary
            
            % Perform a SVD S = V_1 D V_2^T
            [V1_mtx,~,V2_mtx] = svd(S_mtx);
    
            % Extract unitary matrix
            U_mtx = V1_mtx * V2_mtx'; %[U_mtx, ~, ~] = poldecomp(S_mtx); %S_mtx / sqrtm(S_mtx.' * S_mtx); 
            idx = [1:2:2*n, 2:2:2*n];
            U_mtx = U_mtx(idx,idx); % to revert back to normal symplectic matrices

            % Create the complex version
            C_mtx = U_mtx(1:n,1:n) + 1i * U_mtx(n+1:2*n,1:n);
            C_mtx_hist(k, :, :) = C_mtx;
            
        else
            C_mtx = S_mtx;
        end


        % Find eigenvalues
        lambda_vec = eig(C_mtx); lambda_vec = lambda_vec(:).';

        if k > 1
            lambda_vec = match_eigs_to_previous(lambda_vec_hist(k-1,:), lambda_vec);
        end

        lambda_vec_hist(k, :) = lambda_vec;

        theta_vec = imag(log(lambda_vec));
        theta_vec_hist(k, :) = theta_vec;

        % Compute determinant
        d = det(C_mtx);

        % Compute phase [deg]
        if unitary 
            theta_hist(k) = rad2deg(imag(log(d)));
        else
            theta_hist(k) = rad2deg(sum(theta_vec));
        end

        % Consistency checks
        if unitary 
            unit_err_hist(k) = norm(C_mtx' * C_mtx - I_n,'fro');
        else
            unit_err_hist(k) = norm(C_mtx' * C_mtx - I_2n,'fro');
        end
        
        det_err_hist(k) = abs(abs(d) - 1);

    end

    % Unwrap phase history for theta_hist
    theta_hist = rad2deg(unwrap(deg2rad(theta_hist)));

    % Unwrap phase history for thet_vec_hist
    for i = 1 : size(theta_vec_hist,2)
        theta_vec = theta_vec_hist(:,i);
        theta_vec_hist(:,i) = rad2deg(unwrap(theta_vec));
    end

    % Compute total winding angle [deg]
    tot_angle = theta_hist(end) - theta_hist(1);

    % Output
    winding.theta_tot_hist = theta_hist;
    winding.theta_vec_hist = theta_vec_hist;
    winding.lambda_vec_hist = lambda_vec_hist;
    winding.C_mtx_hist = C_mtx_hist;
    winding.theta_tot_end = tot_angle;
    winding.unit_err_hist = unit_err_hist;
    winding.det_err_hist = det_err_hist;

    function lambda_vec_sorted = match_eigs_to_previous(lambda_vec_prev, lambda_vec_now)

        m = length(lambda_vec_prev);
    
        lambda_vec_sorted = zeros(size(lambda_vec_now));
        used_idx = false(size(lambda_vec_now));
    
        for l = 1:m
    
            dist_vec = abs(lambda_vec_now - lambda_vec_prev(l));
            dist_vec(used_idx) = inf;
    
            [~, idx_min] = min(dist_vec);
    
            lambda_vec_sorted(l) = lambda_vec_now(idx_min);
            used_idx(idx_min) = true;
    
        end
    
    end

end



% 
% function theta = compute_winding_number(S_mtx_hist, unitary)
% %==========================================================================
% %
% % Computes the determinant phase winding associated with a time history of
% % symplectic matrices.
% %
% % For each symplectic matrix S(t), the orthogonal/unitary factor is extracted
% % using the polar decomposition through the Singular Value Decomposition
% % (SVD):
% %
% %   S = V_1 D V_2^T
% %
% % so that
% %
% %   U = V_1 V_2^T.
% %
% % Since U is orthogonal-symplectic, it can be identified with the complex
% % unitary matrix
% %
% %   C = U_11 + i U_21.
% %
% % The eigenvalues of C are also computed, providing the individual phase
% % angles associated with each unitary mode:
% %
% %   lambda_i = exp(i theta_i),
% %
% % where theta_i is obtained from the principal value of the complex
% % logarithm.
% %
% % The determinant
% %
% %   d(t) = det(C(t))
% %
% % lies on the unit circle. Its phase equals the sum of the eigenvalue
% % phases, and the winding angle is obtained by unwrapping this determinant
% % phase:
% %
% %   WN = theta(t_f) - theta(t_0).
% %
% % OUTPUT STRUCTURE:
% %
% %   theta:
% %       .hist              - Unwrapped determinant phase history     [deg]
% %       .vec_hist          - Eigenvalue phase history                [deg]
% %       .tot_angle         - Total winding angle                     [deg]
% %       .unit_err_hist     - ||C^*C - I||_F history                  [-]
% %       .det_err_hist      - |||det(C)| - 1| history                 [-]
% %
% % Author: G. Montseny
% % Date: June 18, 2026
% %
% % INPUT:                    Description                          Units
% %
% %   S_mtx_hist     - Symplectic matrix history (N_t x 2n x 2n)   [-]
% %
% % OUTPUT:                   Description                          Units
% %
% %   theta          - Winding angle diagnostic structure          [-]
% %
% %==========================================================================
% 
%     % Initialization
%     N_t = size(S_mtx_hist,1);
%     n = size(S_mtx_hist,2)/2;
% 
%     theta_hist = zeros(N_t,1);
%     theta_vec_hist = zeros(N_t, n);
%     C_mtx_hist = zeros(N_t, n, n);
%     unit_err_hist = zeros(N_t,1);
%     det_err_hist = zeros(N_t,1);
% 
%     I_n = eye(n);
% 
%     % Loop over history
%     for k = 1:N_t
% 
%         % Extract symplectic matrix
%         S_mtx = squeeze(S_mtx_hist(k,:,:));
% 
%         if unitary
%             % Perform a SVD S = V_1 D V_2^T
%             [V1_mtx,~,V2_mtx] = svd(S_mtx);
% 
%             % Extract unitary matrix
%             U_mtx = V1_mtx * V2_mtx'; %[U_mtx, ~, ~] = poldecomp(S_mtx); %S_mtx / sqrtm(S_mtx.' * S_mtx); 
%             idx = [1:2:2*n, 2:2:2*n];
%             U_mtx = U_mtx(idx,idx); % to revert back to normal symplectic matrices
%         else
%             U_mtx = S_mtx;
%         end
% 
%         % Create the complex version
%         C_mtx = U_mtx(1:n,1:n) + 1i * U_mtx(n+1:2*n,1:n);
%         C_mtx_hist(k, :, :) = C_mtx;
% 
%         % Find eigenvalues
%         lambda_vec = eig(C_mtx); lambda_vec = lambda_vec(:).';
%         theta_vec = imag(log(lambda_vec));
%         theta_vec_hist(k, :) = theta_vec;
% 
%         % Compute determinant
%         d = det(C_mtx);
% 
%         % Compute phase [deg]
%         theta_hist(k) = rad2deg(imag(log(d)));
% 
%         % Consistency checks
%         unit_err_hist(k) = norm(C_mtx' * C_mtx - I_n,'fro');
%         det_err_hist(k) = abs(abs(d) - 1);
% 
%     end
% 
%     % Unwrap phase history for theta_hist
%     theta_hist = rad2deg(unwrap(deg2rad(theta_hist)));
% 
%     % Unwrap phase history for thet_vec_hist
%     for i = 1 : n
%         theta_vec = theta_vec_hist(:,i);
%         theta_vec_hist(:,i) = rad2deg(unwrap(theta_vec));
%     end
% 
%     % Compute total winding angle [deg]
%     tot_angle = theta_hist(end) - theta_hist(1);
% 
%     % Output
%     theta.hist = theta_hist;
%     theta.vec_hist = theta_vec_hist;
%     theta.C_mtx_hist = C_mtx_hist;
%     theta.tot_angle = tot_angle;
%     theta.unit_err_hist = unit_err_hist;
%     theta.det_err_hist = det_err_hist;