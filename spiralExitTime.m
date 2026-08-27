function race = spiralExitTime(tSec, elRad, aRad, bRad, hw)
%==========================================================================
% spiralExitTime  Time-marching race between an outward Archimedean spiral
%                 and a time-varying elliptical FOU.
%
%   The gimbal tracks the mean trajectory and, starting at t = tSec(1) from
%   the centre, superimposes an Archimedean spiral that winds OUTWARD. Its
%   radius r_spiral(t) grows at the kinematically/PRF-limited rate. Meanwhile
%   the 3-sigma FOU has semi-axes a(t) (major) and b(t) (minor) that evolve
%   along the pass (they grow toward culmination, then shrink). This routine
%   integrates r_spiral(t) and finds when the spiral crosses the FOU boundary.
%
%   Speed model (quasi-steady, worst-case spiral phase, azimuth keyhole):
%     on-sky spiral speed v(el,r) = min of
%       - PRF no-gap    : v <= ds*PRF
%       - axis velocity : v <= min( Vel , Vaz*cos(el) )   (az stretched 1/cos el)
%       - centripetal   : v <= min( sqrt(Ael*r) , sqrt(Aaz*r*cos(el)) )
%     with ds = 2*Rbeam*(1-alpha) the arm pitch = along-arm pulse spacing.
%   The spiral advances by  dphi = v*dt / sqrt(r^2 + b_c^2),  r = b_c*phi,
%   b_c = ds/2pi.
%
%   Given (all 1xN along the pass, from t = first sample >= start elevation):
%     tSec   sample times [s]
%     elRad  mean-pointing elevation [rad]
%     aRad   FOU semi-major (3-sigma) [rad]
%     bRad   FOU semi-minor (3-sigma) [rad]
%     hw     .Vaz .Vel [deg/s], .Aaz .Ael [deg/s^2], .PRF [Hz],
%            .beamArcsec, .overlap
%
%   Returned struct race:
%     .r            1xN  spiral radius r_spiral(t) [rad]
%     .idxFull .tFull .elFull   first epoch (index/time[s]/elev[deg]) with
%                               r_spiral >= a(t)  (whole ellipse enclosed)
%     .idxFirst .tFirst .elFirst first epoch with r_spiral >= b(t)
%                               (spiral first poking outside the ellipse)
%     .coveredFull  logical, did r_spiral reach a(t) within the pass
%     .dsArcsec     arm pitch [arcsec]
%==========================================================================

    deg = pi/180;  as = deg/3600;
    Vaz = hw.Vaz*deg; Vel = hw.Vel*deg;
    Aaz = hw.Aaz*deg; Ael = hw.Ael*deg;
    if ~isfield(hw,'overlap') || isempty(hw.overlap), hw.overlap = 0; end

    Rbeam = 0.5*hw.beamArcsec*as;
    ds    = 2*Rbeam*(1 - hw.overlap);
    Vsky  = ds*hw.PRF;
    bc    = ds/(2*pi);

    N   = numel(tSec);
    r   = zeros(1, N);
    rr  = Rbeam;                 % start: first footprint already covers centre
    phi = rr/bc;

    for k = 1:N
        r(k) = rr;
        if k == N, break; end
        dt = tSec(k+1) - tSec(k);
        el = elRad(k);  ce = max(cos(el), 1e-3);

        vPRF = Vsky;
        vVel = min(Vel, Vaz*ce);                       % az keyhole (velocity)
        vAcc = min(sqrt(Ael*rr), sqrt(Aaz*rr*ce));     % centripetal (keyhole)
        v    = min([vPRF, vVel, vAcc]);

        dphi = v*dt / sqrt(rr^2 + bc^2);
        phi  = phi + dphi;
        rr   = bc*phi;
    end

    % ---- crossings with the (time-varying) FOU boundary -----------------
    % First crossing: spiral radius first reaches the major axis a(t)
    % (i.e. the swept disk first encloses the whole ellipse).
    kFull  = find(r >= aRad, 1, 'first');
    % First poke outside: spiral radius first reaches the minor axis b(t).
    kFirst = find(r >= bRad, 1, 'first');
    % Persistent coverage: first epoch after which r >= a(t) for ALL later
    % epochs (the FOU never overtakes the spiral again). For a small beam the
    % growing FOU re-engulfs the slow spiral near culmination, so this is much
    % later than the first crossing; for a large beam the two coincide.
    below = r < aRad;
    kLastBelow = find(below, 1, 'last');
    if isempty(kLastBelow)
        kPer = 1;
    elseif kLastBelow < N
        kPer = kLastBelow + 1;
    else
        kPer = [];                         % never persistently covered in pass
    end

    race.r         = r;
    race.dsArcsec  = ds/as;
    race.coveredFull = ~isempty(kFull);
    race.persistCovered = ~isempty(kPer);

    race = setCross(race, 'Full',  kFull,  tSec, elRad, deg);
    race = setCross(race, 'First', kFirst, tSec, elRad, deg);
    race = setCross(race, 'Per',   kPer,   tSec, elRad, deg);
end

% -------------------------------------------------------------------------
function race = setCross(race, tag, k, tSec, elRad, deg)
    if isempty(k) || isnan(k)
        race.(['idx' tag]) = NaN; race.(['t' tag]) = NaN; race.(['el' tag]) = NaN;
    else
        race.(['idx' tag]) = k;
        race.(['t'  tag])  = tSec(k) - tSec(1);
        race.(['el' tag])  = elRad(k)/deg;
    end
end
