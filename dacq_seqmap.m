function C = dacq_seqmap(n)
%==========================================================================
% dacq_seqmap  Perceptually-ordered sequential colormap (viridis anchors).
%   Used to encode TIME, so it must be monotone in lightness: a diverging or
%   hue-cycling map (jet) would read as if the scan doubled back on itself.
%==========================================================================
    A = [0.267 0.005 0.329
         0.229 0.322 0.545
         0.128 0.567 0.551
         0.369 0.789 0.383
         0.993 0.906 0.144];
    x = linspace(0,1,size(A,1));  q = linspace(0,1,n);
    C = [interp1(x,A(:,1),q).', interp1(x,A(:,2),q).', interp1(x,A(:,3),q).'];
end
