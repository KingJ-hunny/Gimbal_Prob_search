function pass = dacq_findpass(cfg)
%==========================================================================
% dacq_findpass  Pick the orbit geometry that yields the requested pass.
%
%   Searches a (RAAN, M0) grid for a circular 800 km orbit so that a pass
%   over the station culminates at ~cfg.elPeakDeg roughly cfg.dtPropHr after
%   epoch. The search uses the SAME RK4+J2 integrator as the sample cloud,
%   so the geometry it reports is the geometry the simulation will see.
%
%   Return pass: .raanDeg .M0Deg .r0 .v0 (state at epoch)
%                .tPeak [s from epoch] .elPeakDeg .tRise .tSet (el = elMask)
%==========================================================================
    mu = 398600.4418;
    aKm = 6378.137 + cfg.altKm;
    e   = cfg.ecc;
    inc = deg2rad(cfg.incDeg);

    raanGrid = deg2rad(0:6:354);
    M0Grid   = deg2rad(0:6:354);
    [RA, MM] = meshgrid(raanGrid, M0Grid);
    RA = RA(:).'; MM = MM(:).';            % 1xK candidates
    K  = numel(RA);

    [r, v] = coe2rvVec(aKm, e, inc, RA, zeros(1,K), MM, mu);

    % ---- coast to the start of the search window -----------------------
    tWin0 = max(900, (cfg.dtPropHr - 1) * 3600);   % search around dtProp, >0
    tWin1 = (cfg.dtPropHr + 1) * 3600;
    hStep = 10;
    [r, v] = dacq_rk4J2(r, v, hStep, round(tWin0/hStep));

    % ---- sweep the window, tracking each candidate's peak elevation ----
    dtW   = 5;
    tS    = tWin0:dtW:tWin1;
    elMax = -ones(1,K); tMax = zeros(1,K);
    for k = 1:numel(tS)
        if k > 1
            [r, v] = dacq_rk4J2(r, v, dtW, 1);
        end
        gm = gmstRad(cfg.epochJD + tS(k)/86400);
        [~, el] = eciToAzElVec(r, gm, cfg.station);
        up = el(:).' > elMax;
        elMax(up) = el(up); tMax(up) = tS(k);
    end

    % ---- pick: peak nearest the requested elevation, inside the window --
    elTgt = deg2rad(cfg.elPeakDeg);
    ok    = tMax > tWin0 + 600 & tMax < tWin1 - 600;    % a real interior peak
    % elevation match dominates; a mild term pulls the pass toward dtPropHr
    cost  = abs(elMax - elTgt) ...
          + 0.020*abs(tMax - cfg.dtPropHr*3600)/3600 + 1e3*(~ok);
    [~, ib] = min(cost);

    pass.raanDeg = rad2deg(RA(ib));
    pass.M0Deg   = rad2deg(MM(ib));
    pass.tPeak   = tMax(ib);
    pass.elPeakDeg = rad2deg(elMax(ib));
    [pass.r0, pass.v0] = coe2rvVec(aKm, e, inc, RA(ib), 0, MM(ib), mu);
    pass.aKm = aKm; pass.ecc = e; pass.incDeg = cfg.incDeg;
end

% =========================================================================
function [r, v] = coe2rvVec(a, e, i, Om, w, M, mu)
% Classical elements -> ECI state, vectorised over 1xK element vectors.
    E = M;
    for it = 1:60
        dE = (E - e*sin(E) - M) ./ (1 - e*cos(E));
        E  = E - dE;
        if max(abs(dE)) < 1e-13, break; end
    end
    nu = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
    p  = a*(1 - e^2);
    rm = p ./ (1 + e*cos(nu));
    xp = rm.*cos(nu);          yp = rm.*sin(nu);
    sp = sqrt(mu/p);
    xd = sp*(-sin(nu));        yd = sp*(e + cos(nu));

    cO=cos(Om); sO=sin(Om); ci=cos(i); si=sin(i); cw=cos(w); sw=sin(w);
    q11 = cO.*cw - sO.*sw*ci;  q12 = -cO.*sw - sO.*cw*ci;
    q21 = sO.*cw + cO.*sw*ci;  q22 = -sO.*sw + cO.*cw*ci;
    q31 = sw*si;               q32 = cw*si;

    r = [q11.*xp + q12.*yp;  q21.*xp + q22.*yp;  q31.*xp + q32.*yp];
    v = [q11.*xd + q12.*yd;  q21.*xd + q22.*yd;  q31.*xd + q32.*yd];
end
