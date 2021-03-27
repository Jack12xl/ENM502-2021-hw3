%% let's test newton method
res = 30;
% U_init = rand(res, res);
m = 2; n = 1;
A = 0.1;
U_init = GuessInit([res, res], A, m, n);
title_str = sprintf('U0 init on %d x %d Grid', res, res);
drawContour(U_init, title_str);

tol = 1e-6;
max_it = 24;
% max_lmbd = 5 * pi^2;
max_lmbd = 60;
min_lmbd = 0;
lmbd_step = 0.1;
ARC_step = 0.1;
ARC_iter_step = 0.2;

lmbd_0 = (m^2 + n^2) * pi^2;
bd_idxes = getBoundaryIdxes([res, res]);
%% solve u0
U_0 = myNewton(res, U_init, lmbd_0, tol, max_it, 2);
U_0 = reshape(U_0, res, res);

title_str = sprintf('U0 on %d x %d Grid, lambda: %d', res, res, lmbd_0);
drawContour(U_0, title_str);
%% Guessing U_1 on lmbd_1 with U_0 and lmbd_0
lmbd_1 = lmbd_0 + lmbd_step;
U_init_lmbd_1 = AnalyticInit(U_0, lmbd_0, lmbd_1, res);
title_str = sprintf('U_init_lmbd_1 on %d x %d Grid, lambda: %d', res, res, lmbd_1);
drawContour(U_init_lmbd_1, title_str);

U_1 = myNewton(res, U_init_lmbd_1, lmbd_1, tol, max_it, 2);

U_1 = reshape(U_1, res, res);

title_str = sprintf('U1 on %d x %d Grid, lambda: %d', res, res, lmbd_1);
drawContour(U_1, title_str);


%% get arc init

[J_hat, b, d_s] = ARCInit(res, U_1, U_0, lmbd_1, lmbd_0);
x = J_hat \ b;
d_u_d_s = x(1:end-1);
d_lmbd_d_s = x(end);

U_2_init = U_1 + ARC_step * reshape(d_u_d_s, res, res);
lmbd_2_init = lmbd_1 + ARC_step * d_lmbd_d_s;

% title_str = sprintf('U2 init on %d x %d Grid, lambda: %d', res, res, lmbd_2_init);
% drawContour(U_2_init, title_str);

[U_arc, lmbd_arc] = myNewton_ARC(res, U_2_init, U_1, lmbd_2_init, lmbd_1, d_s, tol, max_it, 2);
title_str = sprintf('U arc on %d x %d Grid, lambda: %d', res, res, lmbd_arc);
drawContour(U_2_init, title_str);

%% iter
% input
% d_s, U_cur, U_prv, lmbd_cur, lmbd_prv
U_cur = U_arc; U_prv = U_1;
lmbd_cur = lmbd_arc; lmbd_prv = lmbd_1;

iter = 1;
norm_arr = []; lmbd_arr = [];
while (lmbd_cur <= max_lmbd) && (lmbd_cur >= min_lmbd) && (iter <= 32)
    [J, b, ~] = ARCInit(res, U_cur, U_prv, lmbd_cur, lmbd_prv, d_s);
    x = J \ b;
    d_u_d_s = x(1:end-1);
    d_lmbd_d_s = x(end);
    
%     U_init = U_cur + ARC_iter_step * reshape(d_u_d_s, res, res);
%     lmbd_init = lmbd_cur + ARC_iter_step * d_lmbd_d_s;
    U_init = U_cur + d_s * reshape(d_u_d_s, res, res);
    lmbd_init = lmbd_cur + d_s * d_lmbd_d_s;
    
%     title_str = sprintf('U init: %d, lambda: %d', norm(lmbd_init), lmbd_cur);
%     drawContour(U_init, title_str);
    
    % buffer to store previous results
    U_tmp = U_cur; lmbd_tmp = lmbd_cur;
    [U_cur, lmbd_cur] = myNewton_ARC(res, U_init, U_prv, lmbd_init, lmbd_prv, d_s, tol, max_it, 2);
%     U_cur = U_cur + 1e-6;
    U_cur(bd_idxes) = 0;
    % update previous results
    U_prv = U_tmp; lmbd_prv = lmbd_tmp;
    
    if (lmbd_cur <= min_lmbd) || (lmbd_cur >= max_lmbd)
        break;
    end
    
    norm_arr(iter) = norm(U_cur);
    lmbd_arr(iter) = lmbd_cur;
    title_str = sprintf('U norm: %d, lambda: %d', norm_arr(iter), lmbd_cur);
    drawContour(U_cur, title_str);
    
    iter = iter + 1;
end
%%
plot(lmbd_arr, norm_arr, '.');