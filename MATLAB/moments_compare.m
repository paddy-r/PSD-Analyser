% To compare measured and lognormal moments of distribution
function [m_PSD m_lognormal] = moments_compare(x,P,M,S,n,m)
    m_PSD = moment(x,P,n,m);
    m_lognormal = lognormal_moment(M,S,n,m);
end