function[A, b] = NonLinearBVP(U, n, lmbd, ARC_CONT)
%% input
% 
%% output
% @A  J(u)
% @b  -R(u)
    if nargin < 4 || isempty(ARC_CONT)
        ARC_CONT = false;
    end

    h = 1 / n;
    h2 = h.^2;
    h2inv = 1 / h2;
    
    n2 = n.^ 2;
    
    bd_idx = getBoundaryIdxes([n, n]);
    
%% A: the Jacobin matrix A_ij = dR / du
    u = reshape(U, [n2, 1]);
    
    diagArr = - 4.0 * h2inv + lmbd * (1.0 + 2 * u);
    diagArr(bd_idx) = 1;
    
    diagArrSub = ones(n2, 1) * h2inv;
    diagArrSub(bd_idx) = 0;
    
    diagArrLeft = diagArrSub(1+1:end);
    diagArrRight = diagArrSub(1:end-1);
    diagArrTop = diagArrSub(1:end-n);
    diagArrBot = diagArrSub(1+n:end);
    
    A = diag(diagArr) + diag(diagArrLeft, -1) + diag(diagArrRight, 1) + diag(diagArrTop, n) + diag(diagArrBot, -n);
    
    if ARC_CONT
        
    end
%% b: -R(u)
    B2 = U * ( -4 * h2inv + lmbd ) + U.^2 * lmbd;
    
    k = [[0, 1, 0]; [1, 0, 1]; [0, 1, 0]];
    B1 = conv2(U, k, 'same') * h2inv;
    
    b = -(B1 + B2);
    b(bd_idx) = 0;
    
    b = reshape(b, n2, 1);
end