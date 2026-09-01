function [ex, ey] = ellipsePts(Sig, nSig, nPt)
%==========================================================================
% ellipsePts  Boundary points of an nSig-sigma error ellipse of a 2x2
%             covariance Sig, in the same angular units as sqrt(Sig).
%
%   Given:  Sig  2x2 covariance; nSig sigma multiple; nPt point count (opt)
%   Return: ex,ey  1xnPt ellipse boundary coordinates
%==========================================================================
    if nargin < 3, nPt = 72; end
    th = linspace(0, 2*pi, nPt);
    [V, D] = eig(Sig);
    lam = max(diag(D), 0);
    ex = nSig*(V(1,1)*sqrt(lam(1))*cos(th) + V(1,2)*sqrt(lam(2))*sin(th));
    ey = nSig*(V(2,1)*sqrt(lam(1))*cos(th) + V(2,2)*sqrt(lam(2))*sin(th));
end
