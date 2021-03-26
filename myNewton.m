function[U_nxt] = myNewton(n, U_cur, lmbd, tol, max_it, VERBOSE)
    
    if nargin < 4 || isempty(tol)
        tol = 1e-4;
    end
    
    if nargin < 5 || isempty(max_it)
        max_it = 10;
    end
    
    if nargin < 6 || isempty(VERBOSE)
        VERBOSE = 0;
    end
    
    for it = 1:max_it
        [A, b] = NonLinearBVP(n, U_cur, lmbd);
        x = A \ b;
        rsdl = norm(x);
        d_U = reshape(x, n, n);
        U_cur = U_cur + d_U;
        
        if (VERBOSE)
            fprintf('Iteration: %d; Residual: %0.6f\n',it,rsdl);
        end
            
        if rsdl < tol
            break
        end
    end

    U_nxt = U_cur;
        
end