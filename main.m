%% let's test newton method
n = 30;
U_init = rand(n, n);
% U_init = zeros(n, n);
% U_init(10:20, 10:20) = 1;
% bd_idx = getBoundaryIdxes

tol = 1e-6;
max_it = 64;

lmbd_0 = 1;
lmbd_1 = 1;
ARC_CONT = false;
%% solve u0
U_0 = myNewton(U_init, n, lmbd_0, ARC_CONT, tol, max_it);

U_0 = reshape(U_0, n, n);
figure();
contourf(U_0);
colorbar;
title_str = sprintf('U0 on %d x %d Grid, lambda: %d', n, n, lmbd_0);
title(title_str);

% lmbd_1 = 2;
[J_0, ~] = NonLinearBVP(U_0, n, lmbd_0, ARC_CONT);
minus_R_lmbd_0 = - (U_0.^2 + U_0);
minus_R_lmbd_0 = reshape(minus_R_lmbd_0, n^2, 1);
delta_U_lmbd_0 = J_0 \ minus_R_lmbd_0;
delta_U_lmbd_0 = reshape(delta_U_lmbd_0, n, n); 
% get init guess
U_init_lmbd_1 = U_0 + delta_U_lmbd_0 * (lmbd_1 - lmbd_0);
U_1 = myNewton(U_init_lmbd_1, n, lmbd_0, ARC_CONT, tol, max_it);

U_1 = reshape(U_1, n, n);
figure();
contourf(U_1);
colorbar;
title_str = sprintf('U1 on %d x %d Grid, lambda: %d', n, n, lmbd_1);
title(title_str);