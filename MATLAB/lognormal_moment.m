% To get moment(s) of distribution from lognormal parameters of
% distribution
function result = lognormal_moment(varargin)
    M = varargin{1};
    S = varargin{2};
    n = varargin{3};
    if nargin == 4
        m = varargin{4};
    else
        m = 0;
    end

    M1 = exp(n*M) + 0.5*n^2*S^2;
    if m == 0
        M2 = 1.0;
    else
        M2 = exp(m*M) + 0.5*m^2*S^2;
    end
    result = M1/M2;
end