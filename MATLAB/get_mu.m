% To get lognormal mu from M and S parameters
%   mu = exp(M + S^2/2)
function mu = get_mu(M,S)
    mu = exp(M + S^2/2);
end