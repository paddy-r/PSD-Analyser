% Estimate packing fraction of particle species;
%   assumes non-interacting, spherical particles with lognormal sizes
%   - Farr (2013), DOI: https://doi.org/10.1016/j.powtec.2013.04.009
function phi = packing_fraction(M,S)
    phi = 1.0 - 0.57*exp(-S) + 0.2135*exp(-0.57*S/0.2135) + 0.0019*(cos(2.0*pi*(1 - exp(-0.75*S^(0.7) - 0.025*S^4))) - 1);
end