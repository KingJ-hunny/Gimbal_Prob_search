function repIQR(name, v)
%  Report median [Q1, Q3] of a trial vector, NaN-safe.
    w = v(~isnan(v));
    if isempty(w)
        fprintf('     %-14s : n/a (never reached)\n', name);
    else
        fprintf('     %-14s : %8.3f  [%.3f, %.3f]   (n=%d/%d)\n', ...
                name, median(w), qtl(w,0.25), qtl(w,0.75), numel(w), numel(v));
    end
end
