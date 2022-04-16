% To fit to (normalised) lognormal PDF
function [M,S] = fit_lognormal_PDF(x,P,p0)
    B = fminsearch(@(b) norm(P - lognormal_PDF(x,b(1),b(2))),p0);
    M = B(1);
    S = B(2);
end