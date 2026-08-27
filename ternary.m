function v = ternary(cond, a, b)
%TERNARY  Inline conditional: returns a if cond is true, else b.
    if cond, v = a; else, v = b; end
end
