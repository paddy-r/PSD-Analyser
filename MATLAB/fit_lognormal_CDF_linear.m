% To get lognormal parameters M and S CDF
% by linear regression using linearised CDF equation, i.e.
%   ln d  = S.sqrt(2) * inv_err(2C - 1) + M
%   {__}                {_____________}
%    y    =     m              x        + c
% where M = m and
%       S = m / sqrt(2)
function [M,S] = fit_lognormal_CDF_linear(d,C)
    % Transform d,C to X,Y
    X = erfinv(2*C - 1);
    Y = log(d);
    % Linear fit
    p_logn = polyfit(X,Y,1);
    S = p_logn(1)/sqrt(2);
    M = p_logn(2);
end