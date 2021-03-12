%% let's test newton method
n = 30;
U_init = rand(n, n);

tol = 1e-6;
max_it = 64;

lmbd = 10;
ARC_CONT = false;

U = myNewton(U_init, n, lmbd, false, tol, max_it);

U = reshape(U, n, n);
contourf(U);
colorbar;
title_str = sprintf('%d x %d Grid, lambda: %d', n, n, lmbd);
title(title_str);