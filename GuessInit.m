function[U_init] = GuessInit(s, A, m, n)
%% utilize eigen value problem to set initialized value for U
% DESCRIPTIVE TEXT$%^
    
    
    x = linspace(0,1,s(1));
    y = linspace(0,1,s(2));
%     [X, Y] = meshgrid(1:s(1), 1:s(2));
    [X, Y] = meshgrid(x, y);
%     X = X ./ s(1);
%     Y = Y ./ s(2);
    
    U_init = A *sin(m * pi * X) .* sin(n * pi * Y);
end