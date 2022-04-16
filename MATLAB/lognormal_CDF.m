% Get value of lognormal CDF (cumulative distribution function)
% at x with lognormal parameters M and S
%   CDF: C(d)  = 1/2 * [1 + erf((ln x - M)/S.sqrt(2))]
function C = lognormal_CDF(x,M,S)
    C = 0.5*(1+erf((log(x(:))-M)/(S*sqrt(2.0))));
end