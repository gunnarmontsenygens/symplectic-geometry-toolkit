function winding = compute_winding_number(S_mtx_hist, mode, ordering)
%==========================================================================
%
% Computes phase-based winding diagnostics associated with a time history
% of symplectic matrices.
%
% MODES:
%
%   mode = 'direct'
%
%       Uses the eigenvalues of the full real matrix S(t).
%
%       The total phase diagnostic is computed from
%
%           theta_tot(t) = sum_i arg(lambda_i(S(t))).
%
%       No coordinate reordering is applied in this mode because eigenvalues
%       are invariant under permutation similarity.
%
%   mode = 'polar'
%
%       Computes the orthogonal-symplectic factor U(t) from the polar
%       decomposition
%
%           S = R U,
%
%       using the SVD
%
%           S = V1 D V2^T,
%           U = V1 V2^T.
%
%       Then U is written in block form
%
%           U = [ A  -B
%                 B   A ],
%
%       and converted to the complex matrix
%
%           C = A + iB.
%
%       The total phase diagnostic is computed from
%
%           theta_tot(t) = arg(det(C(t))).
%
% ORDERING:
%
%   ordering = 'qp'
%
%       Matrix is already ordered as
%
%           [q1 q2 ... qn p1 p2 ... pn].
%
%   ordering = 'interleaved'
%
%       Matrix is ordered as
%
%           [q1 p1 q2 p2 ... qn pn].
%
%       Reordering is only applied in mode = 'polar', before computing U,
%       R, and C. The original S_mtx_hist is still stored unchanged.
%
% Author: G. Montseny
% Date: June 30, 2026
%
% INPUT:                    Description                          Units
%
%   S_mtx_hist     - Symplectic matrix history                  [-]
%                    (N_t x 2n x 2n)
%
%   mode           - Winding computation mode                   [-]
%
%                    'direct' : eigenvalues of S
%                    'polar'  : orthogonal-symplectic factor U
%
%   ordering       - Coordinate ordering (polar mode only)      [-]
%
%                    'qp'           : [q1 ... qn p1 ... pn]
%                    'interleaved'  : [q1 p1 ... qn pn]
%
% OUTPUT:                   Description                          Units
%
%   winding        - Winding diagnostic structure               [-]
%
%       .S_mtx_hist
%           Original symplectic matrix history                  [-]
%
%       .S_lambda_vec_hist
%           Eigenvalue histories of S                           [-]
%
%       .S_theta_vec_hist
%           Unwrapped eigenvalue phase histories of S           [deg]
%
%       .theta_tot_hist
%           Unwrapped total winding history                     [deg]
%
%       .theta_tot_end
%           Total winding angle                                 [deg]
%
%       .unit_err_hist
%           ||C^*C - I||_F (polar)
%           ||S^TS - I||_F (direct)                             [-]
%
%       .det_err_hist
%           ||det(C)| - 1| (polar)
%           ||det(S)| - 1| (direct)                             [-]
%
%       The following fields are additionally returned when
%       mode = 'polar':
%
%       .U_mtx_hist
%           Orthogonal-symplectic matrix history                [-]
%
%       .U_lambda_vec_hist
%           Eigenvalue histories of U                           [-]
%
%       .U_theta_vec_hist
%           Unwrapped eigenvalue phase histories of U           [deg]
%
%       .R_mtx_hist
%           Positive-definite factor history                    [-]
%
%       .R_lambda_vec_hist
%           Eigenvalue histories of R                           [-]
%
%       .R_theta_vec_hist
%           Unwrapped eigenvalue phase histories of R           [deg]
%
%       .C_mtx_hist
%           Complex unitary matrix history                      [-]
%
%       .C_lambda_vec_hist
%           Eigenvalue histories of C                           [-]
%
%       .C_theta_vec_hist
%           Unwrapped eigenvalue phase histories of C           [deg]
%
%==========================================================================

    % ---------------------------------------------------------------------
    % Initialization
    % ---------------------------------------------------------------------

    N_t = size(S_mtx_hist,1);
    dim = size(S_mtx_hist,2);
    n = dim/2;

    if mod(dim,2) ~= 0
        error('The input matrices must have even dimension.');
    end

    % Interpret mode
    switch lower(mode)
        case 'direct'
            use_polar = false;

        case 'polar'
            use_polar = true;

        otherwise
            error('Unknown mode. Use ''direct'' or ''polar''.');
    end

    % Interpret ordering
    switch lower(ordering)
        case 'qp'
            reorder = false;

        case 'interleaved'
            reorder = true;

        otherwise
            error('Unknown ordering. Use ''qp'' or ''interleaved''.');
    end

    % Reordering index from interleaved to qp ordering
    idx_qp = [1:2:dim, 2:2:dim];

    % Matrix histories
    S_out_mtx_hist = zeros(N_t, dim, dim);

    S_lambda_vec_hist = zeros(N_t, dim);
    S_theta_vec_hist  = zeros(N_t, dim);

    if use_polar
        U_mtx_hist = zeros(N_t, dim, dim);
        R_mtx_hist = zeros(N_t, dim, dim);
        C_mtx_hist = zeros(N_t, n, n);

        U_lambda_vec_hist = zeros(N_t, dim);
        U_theta_vec_hist  = zeros(N_t, dim);

        R_lambda_vec_hist = zeros(N_t, dim);
        R_theta_vec_hist  = zeros(N_t, dim);

        C_lambda_vec_hist = zeros(N_t, n);
        C_theta_vec_hist  = zeros(N_t, n);
    end

    theta_hist = zeros(N_t,1);
    unit_err_hist = zeros(N_t,1);
    det_err_hist = zeros(N_t,1);

    I_n = eye(n);
    I_2n = eye(dim);

    % ---------------------------------------------------------------------
    % Loop over history
    % ---------------------------------------------------------------------

    for k = 1:N_t

        % Extract original matrix
        S_orig_mtx = squeeze(S_mtx_hist(k,:,:));

        % Store original S exactly as input
        S_out_mtx_hist(k,:,:) = S_orig_mtx;

        % -------------------------------------------------------------
        % Direct eigenvalue data of original S
        % -------------------------------------------------------------

        S_lambda_vec = eig(S_orig_mtx);
        S_lambda_vec = S_lambda_vec(:).';

        if k > 1
            S_lambda_vec = match_eigs_to_previous( ...
                S_lambda_vec_hist(k-1,:), S_lambda_vec);
        end

        S_lambda_vec_hist(k,:) = S_lambda_vec;
        S_theta_vec_hist(k,:) = imag(log(S_lambda_vec));

        % -------------------------------------------------------------
        % Polar/unitary mode
        % -------------------------------------------------------------

        if use_polar

            % Work matrix is allowed to be reordered
            S_work_mtx = S_orig_mtx;

            if reorder
                S_work_mtx = S_work_mtx(idx_qp,idx_qp);
            end

            % SVD-based polar decomposition
            [V1_mtx,~,V2_mtx] = svd(S_work_mtx);

            U_mtx = V1_mtx * V2_mtx';

            % Positive-definite factor in the convention S = R U
            R_mtx =   S_work_mtx * U_mtx';

            % Store U and R
            U_mtx_hist(k,:,:) = U_mtx;
            R_mtx_hist(k,:,:) = R_mtx;

            % Eigenvalues of U
            U_lambda_vec = eig(U_mtx);
            U_lambda_vec = U_lambda_vec(:).';

            if k > 1
                U_lambda_vec = match_eigs_to_previous( ...
                    U_lambda_vec_hist(k-1,:), U_lambda_vec);
            end

            U_lambda_vec_hist(k,:) = U_lambda_vec;
            U_theta_vec_hist(k,:) = imag(log(U_lambda_vec));

            % Eigenvalues of R
            R_lambda_vec = eig(R_mtx);
            R_lambda_vec = R_lambda_vec(:).';

            if k > 1
                R_lambda_vec = match_eigs_to_previous( ...
                    R_lambda_vec_hist(k-1,:), R_lambda_vec);
            end

            R_lambda_vec_hist(k,:) = R_lambda_vec;
            R_theta_vec_hist(k,:) = imag(log(R_lambda_vec));

            % Complex representation C = A + iB
            A_mtx = U_mtx(1:n,1:n);
            B_mtx = U_mtx(n+1:dim,1:n);

            C_mtx = A_mtx + 1i*B_mtx;

            C_mtx_hist(k,:,:) = C_mtx;

            % Eigenvalues of C
            C_lambda_vec = eig(C_mtx);
            C_lambda_vec = C_lambda_vec(:).';

            if k > 1
                C_lambda_vec = match_eigs_to_previous( ...
                    C_lambda_vec_hist(k-1,:), C_lambda_vec);
            end

            C_lambda_vec_hist(k,:) = C_lambda_vec;
            C_theta_vec_hist(k,:) = imag(log(C_lambda_vec));

            % Compute determinant
            d = det(C_mtx);

            % Compute phase [deg]
            theta_hist(k) = rad2deg(imag(log(d)));

            % Consistency checks
            unit_err_hist(k) = norm(C_mtx' * C_mtx - I_n,'fro');
            det_err_hist(k) = abs(abs(d) - 1);

        % -------------------------------------------------------------
        % Direct mode
        % -------------------------------------------------------------

        else

            % Direct total phase from eigenvalues of S
            theta_hist(k) = rad2deg(sum(S_theta_vec_hist(k,:)));

            % Diagnostics on S itself
            d = det(S_orig_mtx);

            unit_err_hist(k) = norm(S_orig_mtx' * S_orig_mtx - I_2n,'fro');
            det_err_hist(k) = abs(abs(d) - 1);

        end

    end

    % ---------------------------------------------------------------------
    % Unwrap total phase history
    % ---------------------------------------------------------------------

    theta_hist = rad2deg(unwrap(deg2rad(theta_hist)));

    % ---------------------------------------------------------------------
    % Unwrap S eigenvalue phase histories
    % ---------------------------------------------------------------------

    for i = 1:size(S_theta_vec_hist,2)
        theta_vec = S_theta_vec_hist(:,i);
        S_theta_vec_hist(:,i) = rad2deg(unwrap(theta_vec));
    end

    % ---------------------------------------------------------------------
    % Unwrap polar quantities
    % ---------------------------------------------------------------------

    if use_polar

        for i = 1:size(U_theta_vec_hist,2)
            theta_vec = U_theta_vec_hist(:,i);
            U_theta_vec_hist(:,i) = rad2deg(unwrap(theta_vec));
        end

        for i = 1:size(R_theta_vec_hist,2)
            theta_vec = R_theta_vec_hist(:,i);
            R_theta_vec_hist(:,i) = rad2deg(unwrap(theta_vec));
        end

        for i = 1:size(C_theta_vec_hist,2)
            theta_vec = C_theta_vec_hist(:,i);
            C_theta_vec_hist(:,i) = rad2deg(unwrap(theta_vec));
        end

    end

    % ---------------------------------------------------------------------
    % Compute total winding angle [deg]
    % ---------------------------------------------------------------------

    theta_tot_end = theta_hist(end) - theta_hist(1);

    % ---------------------------------------------------------------------
    % Output
    % ---------------------------------------------------------------------

    winding.S_mtx_hist = S_out_mtx_hist;
    winding.S_lambda_vec_hist = S_lambda_vec_hist;
    winding.S_theta_vec_hist = S_theta_vec_hist;

    if use_polar

        winding.U_mtx_hist = U_mtx_hist;
        winding.U_lambda_vec_hist = U_lambda_vec_hist;
        winding.U_theta_vec_hist = U_theta_vec_hist;

        winding.R_mtx_hist = R_mtx_hist;
        winding.R_lambda_vec_hist = R_lambda_vec_hist;
        winding.R_theta_vec_hist = R_theta_vec_hist;

        winding.C_mtx_hist = C_mtx_hist;
        winding.C_lambda_vec_hist = C_lambda_vec_hist;
        winding.C_theta_vec_hist = C_theta_vec_hist;

    end

    winding.theta_tot_hist = theta_hist;
    winding.theta_tot_end = theta_tot_end;

    winding.unit_err_hist = unit_err_hist;
    winding.det_err_hist = det_err_hist;

    % =====================================================================
    % Nested function: eigenvalue matching
    % =====================================================================

    function lambda_vec_sorted = match_eigs_to_previous(lambda_vec_prev, lambda_vec_now)

        m = length(lambda_vec_prev);

        lambda_vec_sorted = zeros(size(lambda_vec_now));
        used_idx = false(size(lambda_vec_now));

        for l = 1:m

            dist_vec = abs(lambda_vec_now - lambda_vec_prev(l));
            dist_vec(used_idx) = inf;

            [~,idx_min] = min(dist_vec);

            lambda_vec_sorted(l) = lambda_vec_now(idx_min);
            used_idx(idx_min) = true;

        end

    end

end

% function winding = compute_winding_number(S_mtx_hist, unitary, mtx)
% %==========================================================================
% %
% % Computes phase-based winding diagnostics associated with a time history
% % of symplectic matrices.
% %
% % Depending on the value of the input flag `unitary`, the computation is
% % performed using either:
% %
% %   unitary = true
% %
% %       The orthogonal-symplectic factor obtained from the polar
% %       decomposition
% %
% %           S = U P,
% %
% %       where U is extracted via the Singular Value Decomposition (SVD)
% %
% %           S = V_1 D V_2^T,
% %           U = V_1 V_2^T.
% %
% %   unitary = false
% %
% %       The original symplectic matrix S itself.
% %
% % In either case, the matrix is first transformed to the standard
% % (q,p)-ordered symplectic basis and written in block form
% %
% %       [ A  -B ]
% %       [ B  A ].
% %
% % The associated complex matrix is then defined as
% %
% %       C = A + iB.
% %
% % For orthogonal-symplectic matrices this corresponds to the standard
% % identification
% %
% %       U(n) ≅ Sp(2n,R) ∩ O(2n).
% %
% % The eigenvalues of the complex matrix are computed:
% %
% %       lambda_i = rho_i exp(i theta_i),
% %
% % together with their phase angles
% %
% %       theta_i = Im(log(lambda_i)).
% %
% % The determinant
% %
% %       d = det(C)
% %
% % provides a global phase whose argument equals the sum of the individual
% % eigenvalue phases:
% %
% %       arg(d) = sum_i theta_i.
% %
% % The determinant phase history is unwrapped in time and the total winding
% % angle is computed as
% %
% %       theta_tot_end = theta_tot_hist(t_f) - theta_tot_hist(t_0).
% %
% % unitary = true: winding of the orthogonal-symplectic (unitary) factor.
% % unitary = false: evolution of the eigenvalue arguments of the full symplectic matrix S.
% %
% % OUTPUT STRUCTURE:
% %
% %   winding:
% %
% %       .theta_tot_hist
% %           Unwrapped determinant phase history               [deg]
% %
% %       .theta_vec_hist
% %           Unwrapped eigenvalue phase histories              [deg]
% %
% %       .lambda_vec_hist
% %           Eigenvalue histories of C                         [-]
% %
% %       .C_mtx_hist
% %           Complex matrix histories                          [-]
% %
% %       .theta_tot_end
% %           Total winding angle                               [deg]
% %
% %       .unit_err_hist
% %           ||C^* C - I||_F history                           [-]
% %
% %       .det_err_hist
% %           |||det(C)| - 1| history                           [-]
% %
% % Author: G. Montseny
% % Date: June 24, 2026
% %
% % INPUT:                    Description                          Units
% %
% %   S_mtx_hist     - Symplectic matrix history (N_t x 2n x 2n)   [-]
% %   unitary        - Use orthogonal-symplectic factor from the
% %                    polar decomposition if true; otherwise use
% %                    the original symplectic matrix              [-]
% %
% % OUTPUT:                   Description                          Units
% %
% %   winding        - Winding diagnostic structure               [-]
% %
% %==========================================================================
% 
%     % Initialization
%     N_t = size(S_mtx_hist,1);
%     n = size(S_mtx_hist,2)/2;
% 
%     % Different matrices
%     switch mtx
%         case 'Phi'
%             reorder = false;
%         case 'Psi_tilde'
%             reorder = true;
% 
%     end
% 
% 
%     theta_hist = zeros(N_t,1);
%     unit_err_hist = zeros(N_t,1);
%     det_err_hist = zeros(N_t,1);
% 
%     if unitary
%         lambda_vec_hist = zeros(N_t, n);
%         theta_vec_hist = zeros(N_t, n);
%         C_mtx_hist = zeros(N_t, n, n);
%         U_mtx_hist = zeros(N_t, 2*n, 2*n);
%         R_mtx_hist = zeros(N_t, 2*n, 2*n);
%     else
%         lambda_vec_hist = zeros(N_t, 2*n);
%         theta_vec_hist = zeros(N_t, 2*n);
%         C_mtx_hist = zeros(N_t, 2*n, 2*n);
%     end
% 
%     I_n = eye(n);
%     I_2n = eye(2*n);
% 
% 
% 
%     % Loop over history
%     for k = 1 : N_t
% 
%         % Extract symplectic matrix
%         S_mtx = squeeze(S_mtx_hist(k,:,:));
% 
%         if unitary
% 
%             % Perform a SVD S = V_1 D V_2^T
%             [V1_mtx,~,V2_mtx] = svd(S_mtx);
% 
%             % Extract unitary matrix
%             U_mtx = V1_mtx * V2_mtx'; %[U_mtx, ~, ~] = poldecomp(S_mtx); %S_mtx / sqrtm(S_mtx.' * S_mtx); 
% 
%             if reorder 
%                 idx = [1:2:2*n, 2:2:2*n];
%                 U_mtx = U_mtx(idx,idx); % to revert back to normal symplectic matrices
%             end
% 
%             U_mtx_hist(k, :, :) = U_mtx;
% 
%             % Calculate R and append it to its array
%             R_mtx = S_mtx * U_mtx';
%             R_mtx_hist(k, :, :) = R_mtx;
% 
%             % Create the complex version
%             C_mtx = U_mtx(1:n,1:n) + 1i * U_mtx(n+1:2*n,1:n);
%             C_mtx_hist(k, :, :) = C_mtx;
% 
%         else
%             C_mtx = S_mtx;
%         end
% 
% 
%         % Find eigenvalues
%         lambda_vec = eig(C_mtx); lambda_vec = lambda_vec(:).';
% 
%         if k > 1
%             lambda_vec = match_eigs_to_previous(lambda_vec_hist(k-1,:), lambda_vec);
%         end
% 
%         lambda_vec_hist(k, :) = lambda_vec;
% 
%         theta_vec = imag(log(lambda_vec));
%         theta_vec_hist(k, :) = theta_vec;
% 
%         % Compute determinant
%         d = det(C_mtx);
% 
%         % Compute phase [deg]
%         if unitary 
%             theta_hist(k) = rad2deg(imag(log(d)));
%         else
%             theta_hist(k) = rad2deg(sum(theta_vec));
%         end
% 
%         % Consistency checks
%         if unitary 
%             unit_err_hist(k) = norm(C_mtx' * C_mtx - I_n,'fro');
%         else
%             unit_err_hist(k) = norm(C_mtx' * C_mtx - I_2n,'fro');
%         end
% 
%         det_err_hist(k) = abs(abs(d) - 1);
% 
%     end
% 
%     % Unwrap phase history for theta_hist
%     theta_hist = rad2deg(unwrap(deg2rad(theta_hist)));
% 
%     % Unwrap phase history for thet_vec_hist
%     for i = 1 : size(theta_vec_hist,2)
%         theta_vec = theta_vec_hist(:,i);
%         theta_vec_hist(:,i) = rad2deg(unwrap(theta_vec));
%     end
% 
%     % Compute total winding angle [deg]
%     tot_angle = theta_hist(end) - theta_hist(1);
% 
%     % Output
%     winding.theta_tot_hist = theta_hist;
%     winding.theta_vec_hist = theta_vec_hist;
%     winding.lambda_vec_hist = lambda_vec_hist;
%     winding.C_mtx_hist = C_mtx_hist;
% 
%     if unitary
%         winding.U_mtx_hist = U_mtx_hist;
%         winding.R_mtx_hist = R_mtx_hist;
%     end
% 
%     winding.theta_tot_end = tot_angle;
%     winding.unit_err_hist = unit_err_hist;
%     winding.det_err_hist = det_err_hist;
% 
%     function lambda_vec_sorted = match_eigs_to_previous(lambda_vec_prev, lambda_vec_now)
% 
%         m = length(lambda_vec_prev);
% 
%         lambda_vec_sorted = zeros(size(lambda_vec_now));
%         used_idx = false(size(lambda_vec_now));
% 
%         for l = 1:m
% 
%             dist_vec = abs(lambda_vec_now - lambda_vec_prev(l));
%             dist_vec(used_idx) = inf;
% 
%             [~, idx_min] = min(dist_vec);
% 
%             lambda_vec_sorted(l) = lambda_vec_now(idx_min);
%             used_idx(idx_min) = true;
% 
%         end
% 
%     end
% 
% end
% 


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