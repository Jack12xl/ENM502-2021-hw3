function[A_arc, b_arc, d_s] = ARCInit(n, U_n, U_n_1, lmbd_n, lmbd_n_1)
%% input
% 
%% output
% @A  J(u)
% @b  -R(u)
    n2 = n*n;
    [A, ~] = NonLinearBVP(n, U_n, lmbd_n);
%% equation (5)
    d_s = sqrt((lmbd_n - lmbd_n_1).^2 + norm(U_n - U_n_1).^2);
    d_eita_d_lmbd = -2 * (lmbd_n - lmbd_n_1);
    d_eita_d_u = -2 * (U_n - U_n_1);
    d_eita_d_s = 2 * d_s;
    d_R_d_lmbd = U_n .* (1 + U_n);
% insert 
    A_arc = zeros(size(A,1)+1, size(A,2)+1);
    A_arc(1:size(A,1), 1:size(A,2)) = A;
    
    A_arc(end, 1:end-1) = reshape(d_eita_d_u, 1, n2);
    A_arc(end, end) = d_eita_d_lmbd;
    A_arc(1:end-1, end) = reshape(d_R_d_lmbd, n2, 1);
    
    b_arc = zeros(n2 + 1, 1);
    b_arc(end) = -d_eita_d_s;
end