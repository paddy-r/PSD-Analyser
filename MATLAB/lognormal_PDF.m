% Get value of (normalised) lognormal PDF (probability density function)
% at x with lognormal parameters M and S
%   PDF: P(d) = [1/(x.S.sqrt(2.pi)] * exp[-(ln x - M)^2/(2S^2)]
function P = lognormal_PDF(x,M,S)
    P = (1./(S.*sqrt(2*pi)*x(:))).*exp(-((log(x(:))-M).^2)/(2*(S^2)));
end