function harv = harvestScan(sx, sy, tk, pX, pY, ptime, Rbeam)
%==========================================================================
% harvestScan  Which particles are harvested by a spiral scan.
%
%   Scan point k is at offset (sx(k),sy(k)) and is visited at time tk(k).
%   Particle offsets are supplied on a time grid: pX,pY are NxM, sampled at
%   times ptime (1xM). A particle is HARVESTED the first time a scan point
%   lands within Rbeam of the particle's CURRENT position (nearest grid
%   column to tk(k)). For a STATIONARY case pass M=1 (constant offsets) and
%   any tk; every scan point then uses the single column.
%
%   Efficient interval-grouped ring test.
%
%   Given:  sx,sy 1xK [rad]; tk 1xK [s]; pX,pY NxM [rad]; ptime 1xM [s];
%           Rbeam [rad]
%   Return: harv  Nx1 logical
%==========================================================================
    N = size(pX,1); M = size(pX,2); K = numel(sx);
    harv = false(N,1);
    Rb2 = Rbeam^2;
    rho = sqrt(pX.^2 + pY.^2);                 % NxM particle radii

    if M == 1
        mIdx = ones(1,K);
    else
        dtG  = ptime(2) - ptime(1);
        mIdx = min(max(round((tk - ptime(1))/dtG) + 1, 1), M);   % nearest column
    end

    [ms, ord] = sort(mIdx);
    sxs = sx(ord); sys = sy(ord);
    edges = [find(diff(ms) ~= 0), numel(ms)];
    b0 = 1;
    for e = 1:numel(edges)
        b1 = edges(e); m = ms(b0);
        xk = sxs(b0:b1); yk = sys(b0:b1); b0 = b1 + 1;
        rk = sqrt(xk.^2 + yk.^2);
        ring = find(rho(:,m) >= (min(rk)-Rbeam) & rho(:,m) <= (max(rk)+Rbeam) & ~harv);
        if isempty(ring), continue; end
        cx = pX(ring,m); cy = pY(ring,m);
        D2 = (xk(:) - cx(:).').^2 + (yk(:) - cy(:).').^2;    % P x nring
        hit = any(D2 <= Rb2, 1);
        if any(hit), harv(ring(hit)) = true; end
    end
end
