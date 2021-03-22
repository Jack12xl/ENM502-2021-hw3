%% let's test newton method
res = 30;
% U_init = rand(n, n);
m = 2; n = 1;
U_init = GuessInit([res, res], -0.1, m, n);
title_str = sprintf('U0_init on %d x %d Grid', res, res);
drawContour(U_init, title_str);

tol = 1e-6;
max_it = 128;

lmbd_0 = (m^2 + n^2) * pi^2;
lmbd_1 = lmbd_0 + 0.1;
% ARC_CONT = false;
%% solve u0
U_0 = myNewton(res, U_init, lmbd_0, tol, max_it);
U_0 = reshape(U_0, res, res);

title_str = sprintf('U0 on %d x %d Grid, lambda: %d', res, res, lmbd_0);
drawContour(U_0, title_str);
%% Guessing U_1 on lmbd_1 with U_0 and lmbd_0
[J_0, ~] = NonLinearBVP(res, U_0, lmbd_0);
minus_R_lmbd_0 = - (U_0.^2 + U_0);
minus_R_lmbd_0 = reshape(minus_R_lmbd_0, res^2, 1);
delta_U_lmbd_0 = J_0 \ minus_R_lmbd_0;
delta_U_lmbd_0 = reshape(delta_U_lmbd_0, res, res); 
% get init guess
U_init_lmbd_1 = U_0 + delta_U_lmbd_0 * (lmbd_1 - lmbd_0);
U_1 = myNewton(res, U_init_lmbd_1, lmbd_1, tol, max_it);

U_1 = reshape(U_1, res, res);
figure();
contourf(U_1);
colorbar;
title_str = sprintf('U1 on %d x %d Grid, lambda: %d', res, res, lmbd_1);
title(title_str);



%% get arc init
[J_hat, b, d_s] = ARCInit(res, U_1, U_0, lmbd_1, lmbd_0);
x = J_hat \ b;
d_u_d_s = x(1:end-1);
d_lmbd_d_s = x(end);

U_2_init = U_1 + d_s * reshape(d_u_d_s, res, res);
lmbd_2_init = lmbd_1 + d_s * d_lmbd_d_s;
figure();
contourf(U_2_init);
colorbar;
title_str = sprintf('U2\_init on %d x %d Grid, lambda: %d', res, res, lmbd_2_init);
title(title_str);

[U_arc, lmbd_arc] = myNewton_ARC(res, U_2_init, U_1, lmbd_2_init, lmbd_1, d_s, tol, max_it);

