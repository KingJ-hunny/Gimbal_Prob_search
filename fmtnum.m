function s = fmtnum(v)
%FMTNUM  Format a scalar for table output, printing '--' for NaN.
    if isnan(v), s = '  --'; else, s = sprintf('%.1f', v); end
end
