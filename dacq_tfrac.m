function t = dacq_tfrac(res, frac, tmax)
%==========================================================================
% dacq_tfrac  Time at which a given fraction of the sample points has been
%             harvested. Returns NaN if the fraction is never reached.
%==========================================================================
    M = numel(res.tSorted);
    k = max(1, min(M, ceil(frac*M)));
    t = res.tSorted(k);
    if nargin >= 3 && t >= tmax - 1e-12, t = NaN; end
    if isnan(t), t = NaN; end
end
