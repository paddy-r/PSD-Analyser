% Get x corresp. to pth quantile of lognormal distribution
function x_p = lognormal_quantile(M,S,p)
    x_p = exp(M + sqrt(2*S^2)*erfinv(2*p-1.0));
end