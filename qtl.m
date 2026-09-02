function q = qtl(x, p)
%  Simple order-statistic quantile (Octave base has no `quantile`).
    x = x(~isnan(x)); x = sort(x(:));
    if isempty(x), q = NaN; return; end
    k = max(1, min(numel(x), ceil(p*numel(x))));
    q = x(k);
end
