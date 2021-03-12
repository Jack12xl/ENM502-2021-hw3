function[A, b] = NonLinearBVP(U, n, lmbd)
%% input
% 
%% output
% @A  J(u)
% @b  -R(u)

    h = 1 / n;
    h2 = h.^2;
    h2inv = 1 / h2;
    
    n2 = n.^ 2;
    
    bd_idx = getBoundaryIdxes([n, n]);
    
%% J
    u = reshape(U, [n2, 1]);
    
    diagArr = - 4.0 * h2inv + lmbd * (1.0 + 2 * u);
    diagArr(bd_idx) = 0;
    
    diagArrSub = h2inv * ones(n2 , 1);
    diagArrSub(bd_idxes) = 0;
    
    diagArrLeft = diagArrSub(1+1:end);
    diagArrRight = diagArrSub(1:end-1);
    diagArrTop = diagArrSub(1:end-n);
    diagArrBot = diagArrSub(1+n:end);
    
    A = diag(diagArr) + diag(diagArrLeft, -1) + diag(diagArrRight, 1) + diag(diagArrTop, n) + diag(diagArrBot, -n);
end