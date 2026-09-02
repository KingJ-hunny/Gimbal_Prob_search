function q = quantileish(x, p)
  s = sort(x(:)); k = max(1,min(numel(s),ceil(p*numel(s)))); q = s(k);
end
