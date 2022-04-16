% To get lognormal sigma from M and S parameters
%   sigma = sqrt(exp(S^2 + 2M) * (exp(S^2) - 1))
function sigma = get_sigma(M,S)
    sigma = sqrt(exp(S^2 + 2*M) * (exp(S^2) - 1));
end