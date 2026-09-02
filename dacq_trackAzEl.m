function [AZ, EL, Rng] = dacq_trackAzEl(r0, v0, jd0, tvec, station, hMax)
%==========================================================================
% dacq_trackAzEl  Topocentric az/el history of many states over a time grid.
%
%   Propagates the 3xM state set with RK4+J2 (dacq_rk4J2) and converts to
%   topocentric az/el at every epoch of tvec.
%
%   Given:
%     r0,v0    3xM  states at jd0 [km],[km/s]
%     jd0      scalar Julian Date of the initial state
%     tvec     1xK  epochs, seconds from jd0 (increasing, uniform)
%     station  struct .latDeg .lonDeg .rEcef
%     hMax     max RK4 sub-step [s]
%
%   Return:
%     AZ, EL   MxK  azimuth / elevation [rad]
%     Rng      MxK  slant range [km]
%==========================================================================
    M = size(r0,2); K = numel(tvec);
    AZ = zeros(M,K); EL = zeros(M,K); Rng = zeros(M,K);

    r = r0; v = v0; tPrev = 0;
    for k = 1:K
        dt = tvec(k) - tPrev;
        if dt > 0
            ns = max(1, ceil(dt / hMax));
            [r, v] = dacq_rk4J2(r, v, dt/ns, ns);
            tPrev = tvec(k);
        end
        gm = gmstRad(jd0 + tvec(k)/86400);
        [az, el, rg] = eciToAzElVec(r, gm, station);
        AZ(:,k) = az(:); EL(:,k) = el(:); Rng(:,k) = rg(:);
    end
end
