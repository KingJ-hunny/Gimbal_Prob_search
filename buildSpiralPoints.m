function [sx, sy, K, Rs] = buildSpiralPoints(aMax, Rbeam, alpha)
%==========================================================================
% buildSpiralPoints  Constant-linear-velocity Archimedean spiral scan points
%                    out to radius aMax (the FOU major axis).
%
%       R_s = 2*alpha*Rbeam        (alpha=1 -> adjacent footprints tangent)
%       theta_k = sqrt(4*pi*k),  r_k = (R_s/2pi)*theta_k = R_s*sqrt(k/pi)
%       point k = (r_k cos theta_k, r_k sin theta_k),  k = 1..K
%   K is the number of points until r_K reaches aMax.
%
%   Given:  aMax  scan radius (FOU semi-major) [rad]
%           Rbeam beam footprint radius [rad]
%           alpha overlap factor
%   Return: sx, sy  1xK scan-point offsets [rad]
%           K       point count
%           Rs      point spacing / arm pitch [rad]
%==========================================================================
    Rs = 2*alpha*Rbeam;
    K  = max(1, floor(pi*(aMax/Rs)^2));
    k  = 1:K;
    th = sqrt(4*pi*k);
    r  = (Rs/(2*pi))*th;
    sx = r.*cos(th);
    sy = r.*sin(th);
end
