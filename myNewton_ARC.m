function[U_nxt] = myNewton_ARC(n, U_cur, U_prv, lmbd, lmbd_prv, tol, max_it)
    
    if nargin < 4 || isempty(tol)
        tol = 1e-4;
    end
    
    if nargin < 5 || isempty(max_it)
        max_it = 10;
    end
    
    
    for it = 1:max_it
        [A, b, d_s] = NonLinearBVP_ARC(n, U_cur, lmbd);
        x = A \ b;
        
%         d_u_d_s = x(1:end-1);
%         d_lmbd_d_s = x(end);
        
        
%         U_cur = U_cur + d_s * d_
        
        rsdl = norm(x);
        d_U = reshape(x, n, n);
        U_cur = U_cur + d_U;

        fprintf('Iteration: %d; Residual: %0.6f\n',it,rsdl);

        if rsdl < tol
            break
        end
    end

    U_nxt = U_cur;
        
end