%% let's test newton method
n = 30;
U_init = rand(n, n);

tol = 1e-6;
max_it = 64;

lmbd_0 = 10;
ARC_CONT = false;
%% solve u0
U_0 = myNewton(U_init, n, lmbd_0, false, tol, max_it);

U_0 = reshape(U_0, n, n);
contourf(U_0);
colorbar;
title_str = sprintf('U0 on %d x %d Grid, lambda: %d', n, n, lmbd_0);
title(title_str);

