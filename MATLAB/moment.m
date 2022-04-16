% To get moment(s) of distribution from measured (normalised) PDF
function result = moment(varargin)
    x = varargin{1};
    f = varargin{2};
    n = varargin{3};
    if nargin == 4
        m = varargin{4};
    else
        m = 0;
    end

    M1 = sum((x.^n).*f);
    if m == 0
        M2 = 1.0;
    else
        M2 = sum((x.^m).*f);
    end
    result = M1/M2;
end