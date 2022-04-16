% To fit to lognormal CDF
function [M,S] = fit_lognormal_CDF(x,C,p0)
    B = fminsearch(@(b) norm(C - lognormal_CDF(x,b(1),b(2))),p0);
    M = B(1);
    S = B(2);
end