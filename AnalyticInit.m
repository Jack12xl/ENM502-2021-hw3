function [U_init_lmbd_1] = AnalyticInit(U_0, lmbd_0, lmbd_1, res)
%% Use analytic continuation(first order) to guess init solution at lambda_1 from lambda_0
    [J_0, ~] = NonLinearBVP(res, U_0, lmbd_0);
    minus_R_lmbd_0 = - (U_0.^2 + U_0);
    minus_R_lmbd_0 = reshape(minus_R_lmbd_0, res^2, 1);
    delta_U_lmbd_0 = J_0 \ minus_R_lmbd_0;
    delta_U_lmbd_0 = reshape(delta_U_lmbd_0, res, res); 
    % get init guess
    U_init_lmbd_1 = U_0 + delta_U_lmbd_0 * (lmbd_1 - lmbd_0);
end