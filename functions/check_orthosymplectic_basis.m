function check = check_orthosymplectic_basis(u_hat_set, v_hat_set, w, tol)
%==========================================================================
%
% Checks whether a set of vectors forms an orthosymplectic basis in
% R^{2n}. The function verifies:
%
%   - unit norm
%   - linear independence
%   - Euclidean orthogonality
%   - symplectic orthogonality
%
% with respect to a specified symplectic 2-form.
%
% Author: G. Montseny
% Date: May 20, 2026
%
% INPUTS:                   Description                          Units
%
%   u_hat_set    - Set of u basis vectors in R^{2n}             [-]
%                  [n x 2n]
%
%   v_hat_set    - Set of v basis vectors in R^{2n}             [-]
%                  [n x 2n]
%
%   w            - Symplectic 2-form function handle            [-]
%
%   tol          - Numerical tolerance                          [-]
%
% OUTPUTS:                  Description                          Units
%
%   check        - Structure containing orthosymplectic         [-]
%                  consistency checks
%
%==========================================================================

    % INITIALIZATION
    n = size(u_hat_set, 1) ;
    N = 2 * n ;
    uv_hat_set = [u_hat_set; v_hat_set];
    check.number_basis_vecs = N;

    if nargin < 3 || isempty(w)
        w = @standard_symplectic_2form;
    end

    if nargin < 4 || isempty(tol)
        tol = 1e-10;
    end

    % (0) CHECK IF THE VECTORS ARE UNIT VECTORS
    unit = true;
    for i = 1 : N
        if abs(norm(uv_hat_set(i,:))-1) > tol
            unit = false;
        end
    end

    if unit == true
        check.unit = 'true';
    else
        check.unit = 'false';
    end    

    
    % (I) CHECK IF THE VECTORS ARE LINEARLY INDEPENDENT

    if check.number_basis_vecs == rank(uv_hat_set)
        check.LI = 'true';
    else
        check.LI = 'false';
    end

    % (II) CHECK IF THE VECTORS ARE ORTHOGONAL
    orth = true;

    for i = 1 : N
        for j = 1 : N
            if i ~= j
                if abs(dot(uv_hat_set(i, :), uv_hat_set(j, :))) > tol 
                    orth = false;
                end
            end
        end
    end

    
    if orth == true
        check.orthogonal = 'true';
    else
        check.orthogonal = 'false';
    end


    %(III) CHECK IF THE BASIS IS SYMPLECTIC
    sympl = true;

    for i = 1 : n
        for j = 1 : n

            w_u_ij = w(u_hat_set(i, :), u_hat_set(j, :));
            w_v_ij = w(v_hat_set(i, :), v_hat_set(j, :));
            w_uv_ij = w(u_hat_set(i, :), v_hat_set(j, :));

            if abs(w_u_ij) > tol 
                sympl = false;
            end

            if abs(w_v_ij) > tol 
                sympl = false;
            end

            if i == j
                if abs(w_uv_ij-1) > tol 
                    sympl = false;
                end
            else
                if abs(w_uv_ij) > tol 
                    sympl = false;
                end                
            end

        end
    end

    if sympl == true
        check.symplectic = 'true';
    else
        check.symplectic = 'false';
    end

end