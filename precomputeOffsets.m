function [dX, dY, elMeanDeg] = precomputeOffsets(tvec, tle, station, sigmaRIC, dR)
%==========================================================================
% precomputeOffsets  Co-moving sky offsets of N candidate particles over time.
%
%   Each particle is a fixed RIC offset dR(i,:) [km]. At each time it is placed
%   at r_mean(t)+Q(t)*dR' (ECI), projected to topocentric az/el, and expressed
%   as an offset from the mean pointing in the local sky-tangent frame:
%       dX = cross-elevation offset, dY = elevation offset  [rad]
%   (i.e. the spiral is assumed centred on the mean orbit, so these offsets are
%   what the scan pattern must cover.)
%
%   Given:
%     tvec      1xM   times [s from TLE epoch]
%     tle       parsed TLE struct
%     station   .latDeg .lonDeg .rEcef
%     sigmaRIC  struct (only the RIC basis is used; magnitudes irrelevant here)
%     dR        Nx3   RIC offsets [km]
%
%   Returned:
%     dX, dY    NxM single  local offsets [rad]
%     elMeanDeg 1xM         mean-pointing elevation [deg]
%==========================================================================
    cfg.mode = 'tle'; cfg.tle = tle; cfg.tSec = tvec;
    trk = generateTrack(cfg);
    cov = generateCovariance(trk, sigmaRIC);   % Rhat/Ihat/Chat are geometric

    N = size(dR,1); M = numel(tvec);
    dX = zeros(N,M,'single'); dY = zeros(N,M,'single'); elMeanDeg = zeros(1,M);
    for j = 1:M
        gm = gmstRad(trk.jd(j));
        [az0, el0] = computeAzAlt(trk.rEci(:,j), gm, station);
        Q = [cov.Rhat(:,j), cov.Ihat(:,j), cov.Chat(:,j)];
        candEci = trk.rEci(:,j) + Q*dR.';
        [caz, cel] = eciToAzElVec(candEci, gm, station);
        dA = mod(caz - az0 + pi, 2*pi) - pi;
        dX(:,j) = single((dA .* cos(el0)).');
        dY(:,j) = single((cel - el0).');
        elMeanDeg(j) = rad2deg(el0);
    end
end
