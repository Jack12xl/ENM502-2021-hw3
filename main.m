%% let's test newton method
n = 30;
U_init = rand(n, n);
% U_init = zeros(n, n);
% U_init(10:20, 10:20) = 1;
% bd_idx = getBoundaryIdxes

tol = 1e-6;
max_it = 64;

lmbd_0 = 1;
lmbd_1 = 2;
ARC_CONT = false;
%% solve u0
U_0 = myNewton(n, U_init, lmbd_0, tol, max_it);

U_0 = reshape(U_0, n, n);
figure();
contourf(U_0);
colorbar;
title_str = sprintf('U0 on %d x %d Grid, lambda: %d', n, n, lmbd_0);
title(title_str);

% lmbd_1 = 2;
[J_0, ~] = NonLinearBVP(n, U_0,  lmbd_0);
minus_R_lmbd_0 = - (U_0.^2 + U_0);
minus_R_lmbd_0 = reshape(minus_R_lmbd_0, n^2, 1);
delta_U_lmbd_0 = J_0 \ minus_R_lmbd_0;
delta_U_lmbd_0 = reshape(delta_U_lmbd_0, n, n); 
% get init guess
U_init_lmbd_1 = U_0 + delta_U_lmbd_0 * (lmbd_1 - lmbd_0);
U_1 = myNewton(n, U_init_lmbd_1, lmbd_0, tol, max_it);

U_1 = reshape(U_1, n, n);
figure();
contourf(U_1);
colorbar;
title_str = sprintf('U1 on %d x %d Grid, lambda: %d', n, n, lmbd_1);
title(title_str);

% Assume s at p_0 equal 0
%% equation (5)
% d_s_1 = sqrt((lmbd_1 - lmbd_0).^2 + norm(U_1 - U_0).^2);
% d_eita_d_lmbd = -2 * (lmbd_1 - lmbd_0);
% d_eita_d_u = -2 * (U_1 - U_0);
% d_eita_d_s = 2 * d_s_1;

%% get arc init
[J_hat, b, d_s] = ARCInit(n, U_1, U_0, lmbd_1, lmbd_0);
x = J_hat \ b;
d_u_d_s = x(1:end-1);
d_lmbd_d_s = x(end);

U_2_init = U_1 + d_s * reshape(d_u_d_s, n, n);
lmbd_2_init = lmbd_1 + d_s * d_lmbd_d_s;



