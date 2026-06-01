function u_tilde_vec_set = generate_LI_vec(u1_tilde_vec)
%==========================================================================
%
% Generates n linearly independent vectors in R^{2n} from a prescribed
% nonzero vector. The first vector of the set is the prescribed vector,
% while the remaining n-1 vectors are constructed from standard basis
% vectors chosen such that the resulting set remains linearly independent.
%
% The construction assumes:
%
%   - u1_tilde_vec ~= 0
%   - N = length(u1_tilde_vec) = 2n
%
% Author: G. Montseny
% Date: May 20, 2026
%
% INPUTS:                   Description                          Units
%
%   u1_tilde_vec  - Prescribed nonzero vector                   [-]
%                   in R^{2n} [2n x 1]
%
% OUTPUTS:                  Description                          Units
%
%   u_tilde_vec_set  - Matrix whose rows form a linearly           [-]
%                   independent set of n vectors
%                   in R^{2n} [n x 2n]
%
%==========================================================================


    % Initialization
    u1_tilde_vec = u1_tilde_vec(:);
    N = length(u1_tilde_vec);
    n = N/2;
    u_tilde_vec_set = zeros(n, N);

    % The first u_tilde_vec is trivial
    u_tilde_vec_set(1, :) = u1_tilde_vec';

    % If the first component of the given vector is not zero
    if u1_tilde_vec(1) ~= 0
        
        % Build the remaining u_tilde_vecs
        for i = 2 : n
            ui_tilde_vec = [zeros(1,i-1), 1, zeros(1,N-i)];
            u_tilde_vec_set(i,:) = ui_tilde_vec;
        end
    
    % If the first component is zero
    else
        
        % Introduce new index
        j = 2;
        
        % Find the first component of the given vector that is not zero
        while u1_tilde_vec(j) == 0
            j = j + 1;
        end

        % Build the remaining u_tilde_vecs
        
        u_tilde_vec_set(2,:) = [1, zeros(1,N-1)];

        k = 2;

        for i = 3 : n
            
            if k ~= j
                ui_tilde_vec = [zeros(1,k-1), 1, zeros(1,N-k)];
                k = k + 1;
            else
                k = k + 1;
                ui_tilde_vec = [zeros(1,k-1), 1, zeros(1,N-k)];
            end

            u_tilde_vec_set(i,:) = ui_tilde_vec;
        end

    end


end

