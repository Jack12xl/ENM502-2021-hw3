%% let's test newton method
n = 8;
U_init = rand(8, 8);

tol = 1e-5;
max_it = 16;

lmbd = 4;
ARC_CONT = false;

U = myNewton(U_init, n, lmbd, false, tol, max_it);