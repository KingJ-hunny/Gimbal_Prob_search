function trk = generateTrack(cfg)
%==========================================================================
% generateTrack  (FUNCTION 1)  Mean-orbit track generator.
%
%   Produces the nominal (mean) satellite state history in ECI for a range
%   of orbit regimes.  Two propagation back-ends are supported and selected
%   by cfg.mode:
%
%     'tle'    : realistic mode. SGP4 propagation of a NORAD TLE.
%                Use cfg.line1/cfg.line2 (raw TLE) or cfg.tle (parseTLE out).
%
%     'kepler' : analytic mode for arbitrary regimes (LEO / MEO / GEO /
%                Molniya-type high-eccentricity). Two-body propagation with
%                optional J2 secular drift of RAAN, argument of perigee and
%                mean anomaly.  Use cfg.elements (classical elements).
%
%   -------------------------------------------------------------------
%   Common time specification (choose ONE):
%     cfg.tSec         1xN explicit sample times, seconds from epoch, OR
%     cfg.tSpanSec = [t0 t1] and cfg.dtSec = step   -> t0:dt:t1
%
%   -------------------------------------------------------------------
%   cfg fields for mode 'tle':
%     cfg.line1, cfg.line2   raw TLE strings              (or)
%     cfg.tle                struct from parseTLE
%
%   cfg fields for mode 'kepler':
%     cfg.elements.aKm        semi-major axis        [km]
%     cfg.elements.ecc        eccentricity           [-]
%     cfg.elements.incDeg     inclination            [deg]
%     cfg.elements.raanDeg    RAAN                   [deg]
%     cfg.elements.argpDeg    argument of perigee    [deg]
%     cfg.elements.M0Deg      mean anomaly at epoch  [deg]   (or nu0Deg)
%     cfg.elements.epochJD    epoch Julian Date      [days]  (default J2000)
%     cfg.useJ2               logical, add J2 secular rates  (default true)
%
%   -------------------------------------------------------------------
%   Returned struct trk:
%     .mode      propagation mode used
%     .label     human-readable orbit label (regime + id)
%     .tSec      1xN  sample times, seconds from epoch
%     .jd        1xN  absolute Julian Date of each sample (UTC)
%     .rEci      3xN  ECI position [km]
%     .vEci      3xN  ECI velocity [km/s]
%     .epochJD   scalar epoch Julian Date
%     .periodMin orbital period estimate [min]
%     .tle       (mode 'tle' only) the parsed TLE struct
%==========================================================================

    mu = 398600.4418;                      % Earth GM [km^3/s^2]

    % ---- resolve the time vector -------------------------------------
    if isfield(cfg, 'tSec') && ~isempty(cfg.tSec)
        tSec = cfg.tSec(:).';
    elseif isfield(cfg, 'tSpanSec') && isfield(cfg, 'dtSec')
        tSec = cfg.tSpanSec(1):cfg.dtSec:cfg.tSpanSec(2);
    else
        error('generateTrack:time', 'Provide cfg.tSec OR cfg.tSpanSec+cfg.dtSec.');
    end
    N = numel(tSec);

    rEci = zeros(3, N);
    vEci = zeros(3, N);

    switch lower(cfg.mode)
    % ==================================================================
    case 'tle'
        if isfield(cfg, 'tle') && ~isempty(cfg.tle)
            Tle = cfg.tle;
        else
            Tle = parseTLE(cfg.line1, cfg.line2);
        end
        epochJD = Tle.epochJD;
        for k = 1:N
            [r, v] = sgp4(Tle, tSec(k)/60);        % SGP4 wants minutes-from-epoch
            rEci(:,k) = r(:);
            vEci(:,k) = v(:);
        end
        aApprox   = (mu / (Tle.meanMotion*2*pi/86400)^2)^(1/3);
        periodMin = 1440 / Tle.meanMotion;
        trk.tle   = Tle;
        trk.label = sprintf('%s | SGP4 sat %d', regimeLabel(aApprox, Tle.eccentricity), Tle.satnum);

    % ==================================================================
    case 'kepler'
        el = cfg.elements;
        if ~isfield(el, 'epochJD') || isempty(el.epochJD), el.epochJD = 2451545.0; end
        if ~isfield(cfg, 'useJ2'),  cfg.useJ2 = true; end
        epochJD = el.epochJD;

        a = el.aKm; e = el.ecc;
        i    = deg2rad(el.incDeg);
        raan = deg2rad(el.raanDeg);
        argp = deg2rad(el.argpDeg);
        if isfield(el, 'M0Deg') && ~isempty(el.M0Deg)
            M0 = deg2rad(el.M0Deg);
        else
            M0 = trueToMeanAnomaly(deg2rad(el.nu0Deg), e);
        end

        n = sqrt(mu / a^3);                    % mean motion [rad/s]

        % J2 secular rates (rad/s)
        raanDot = 0; argpDot = 0; M_J2Dot = 0;
        if cfg.useJ2
            J2 = 1.08262668e-3; Re = 6378.137;
            p  = a*(1 - e^2);
            fac = n * J2 * (Re/p)^2;
            raanDot = -1.5 * fac * cos(i);
            argpDot =  0.75 * fac * (5*cos(i)^2 - 1);         % = 1.5*fac*(2 - 2.5 sin^2 i)
            M_J2Dot =  0.75 * fac * sqrt(1 - e^2) * (3*cos(i)^2 - 1);
        end

        for k = 1:N
            t = tSec(k);
            M  = M0   + (n + M_J2Dot)*t;
            Om = raan + raanDot*t;
            w  = argp + argpDot*t;
            [rp, vp] = coe2rv(a, e, i, Om, w, M, mu);
            rEci(:,k) = rp;
            vEci(:,k) = vp;
        end
        periodMin = 2*pi/n / 60;
        trk.label = sprintf('%s | Kepler a=%.0fkm e=%.3f i=%.1fdeg', ...
                            regimeLabel(a, e), a, e, el.incDeg);

    otherwise
        error('generateTrack:mode', 'Unknown cfg.mode "%s" (use ''tle'' or ''kepler'').', cfg.mode);
    end

    trk.mode      = lower(cfg.mode);
    trk.tSec      = tSec;
    trk.jd        = epochJD + tSec/86400;
    trk.rEci      = rEci;
    trk.vEci      = vEci;
    trk.epochJD   = epochJD;
    trk.periodMin = periodMin;
end

% =========================================================================
function [r, v] = coe2rv(a, e, i, Om, w, M, mu)
% Classical elements (with mean anomaly M) -> ECI position/velocity.
    E = keplerE(M, e);
    nu = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
    p  = a*(1 - e^2);
    r_pf_mag = p / (1 + e*cos(nu));

    r_pf = r_pf_mag * [cos(nu); sin(nu); 0];
    v_pf = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];

    cO=cos(Om); sO=sin(Om); ci=cos(i); si=sin(i); cw=cos(w); sw=sin(w);
    Q = [ cO*cw - sO*sw*ci,  -cO*sw - sO*cw*ci,   sO*si;
          sO*cw + cO*sw*ci,  -sO*sw + cO*cw*ci,  -cO*si;
          sw*si,              cw*si,              ci ];
    r = Q * r_pf;
    v = Q * v_pf;
end

% =========================================================================
function E = keplerE(M, e)
% Solve Kepler's equation M = E - e sin E by Newton iteration.
    M = mod(M, 2*pi);
    E = M; if e > 0.8, E = pi; end
    for it = 1:100
        dE = (E - e*sin(E) - M) / (1 - e*cos(E));
        E = E - dE;
        if abs(dE) < 1e-12, break; end
    end
end

% =========================================================================
function M = trueToMeanAnomaly(nu, e)
    E = 2*atan2( sqrt(1-e)*sin(nu/2), sqrt(1+e)*cos(nu/2) );
    M = E - e*sin(E);
end

% =========================================================================
function lab = regimeLabel(aKm, e)
% Coarse orbit-regime label from semi-major axis and eccentricity.
    rp = aKm*(1 - e);  ra = aKm*(1 + e);
    hp = rp - 6378.137; ha = ra - 6378.137;
    if e >= 0.25
        lab = 'HEO';
    elseif ha < 2000
        lab = 'LEO';
    elseif ha < 35000
        lab = 'MEO';
    else
        lab = 'GEO';
    end
end
