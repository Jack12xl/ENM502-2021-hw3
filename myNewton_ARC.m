function[U_nxt, lmbd_nxt] = myNewton_ARC(n, U_cur, U_prv, lmbd, lmbd_prv, d_s, tol, max_it, VERBOSE)
    
    if nargin < 7 || isempty(tol)
        tol = 1e-4;
    end
    
    if nargin < 8 || isempty(max_it)
        max_it = 10;
    end
    
    if nargin < 9 || isempty(VERBOSE)
        % 2: printf iteration norm
        % 1: print result norm
        VERBOSE = 0;
        
    end
    
    
    bd_idxes = getBoundaryIdxes(size(U_cur));
    U_cur(bd_idxes) = 0;
    U_prv(bd_idxes) = 0;
    
    for it = 1:max_it
        [J_hat, minus_R_hat] = NonLinearBVP_ARC(n, U_cur, U_prv, lmbd, lmbd_prv, d_s);
        x = J_hat \ minus_R_hat;
        
        d_U = x(1:end-1);
        d_lmbd = x(end);
        
        lmbd_prv = lmbd;
        U_prv = U_cur;
        
        rsdl = norm(x);
        d_U = reshape(d_U, n, n);
        U_cur = U_cur + d_U;
        lmbd = lmbd + d_lmbd;
        
        U_cur(bd_idxes) = 0;
        
        if (VERBOSE >= 1)
            fprintf('Iteration: %d; Residual: %0.6f\n',it,rsdl);
        end
            
        if rsdl < tol
            break
        end
    end

    U_nxt = U_cur;
    lmbd_nxt = lmbd;
    if (VERBOSE >= 2)
        fprintf('Take %d iterations; Residual: %0.6f\n',it, rsdl);
    end
end