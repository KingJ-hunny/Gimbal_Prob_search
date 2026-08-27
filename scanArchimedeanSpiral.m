function out = scanArchimedeanSpiral(aRad, el0Rad, hw)
%==========================================================================
% scanArchimedeanSpiral  Conventional Archimedean-spiral scan of a circular
%                        FOU, and the time to cover it under real gimbal
%                        velocity/acceleration limits, PRF, and beam width.
%
%   The conventional baseline (Kim et al. "ASS") centres a circular
%   Archimedean spiral on the mean-trajectory pointing and scans outward to
%   r = a (the FOU semi-major axis, so the whole elliptical FOU is enclosed).
%   Arm pitch and along-arm pulse spacing are both
%        ds = 2*Rbeam*(1-alpha) = beamwidth*(1-alpha).
%
%   Time to complete the scan is computed as a time-optimal traversal of the
%   fixed spiral PATH subject to three independent limits:
%     (1) PRF no-gap    : on-sky speed <= ds*PRF  (consecutive footprints
%                         must overlap along the arm, spacing <= ds),
%     (2) gimbal rate   : |d/dt Az| <= Vaz , |d/dt El| <= Vel,
%     (3) gimbal accel  : |d2/dt2 Az| <= Aaz, |d2/dt2 El| <= Ael.
%   Near zenith the azimuth axis is stretched by 1/cos(el0) (keyhole), which
%   is included exactly in the (Az,El) path.
%
%   Given:
%     aRad    FOU semi-major axis (scan-circle radius), on-sky angle [rad]
%     el0Rad  mean-pointing elevation at the scan epoch [rad] (az stretch)
%     hw      struct:
%               .Vaz .Vel   max axis rate  [deg/s]   (e.g. 20, 10)
%               .Aaz .Ael   max axis accel [deg/s^2] (e.g. 5, 2)
%               .PRF        pulse rate     [Hz]      (e.g. 2000)
%               .beamArcsec beam width (footprint diameter) [arcsec]
%               .overlap    alpha in [0,1) (default 0)
%
%   Returned struct out:
%     .T_full     scan time with velocity+accel+PRF limits [s]  (primary)
%     .T_vel      scan time with velocity+PRF only (accel ignored) [s]
%     .T_dwellLB  idealised point-based lower bound K/PRF [s]
%     .K          Archimedean scan-point count pi*a^2/ds^2 + 1
%     .Nturn      number of spiral turns a/ds
%     .LskyDeg    on-sky spiral path length [deg]
%     .dsArcsec   arm pitch / pulse spacing [arcsec]
%     .VskyDeg    PRF-limited on-sky speed [deg/s]
%     .bindFrac   struct, fraction of the path limited by each constraint
%     .phi .phidot .MVC .az .el   path + speed profile (for plotting)
%==========================================================================

    if ~isfield(hw,'overlap') || isempty(hw.overlap), hw.overlap = 0; end

    deg = pi/180;  as = deg/3600;                 % rad per deg, rad per arcsec
    Vaz = hw.Vaz*deg; Vel = hw.Vel*deg;           % rad/s
    Aaz = hw.Aaz*deg; Ael = hw.Ael*deg;           % rad/s^2

    Rbeam = 0.5*hw.beamArcsec*as;                 % footprint radius [rad]
    ds    = 2*Rbeam*(1 - hw.overlap);             % arm pitch & pulse spacing [rad]
    Vsky  = ds*hw.PRF;                            % PRF no-gap on-sky speed cap [rad/s]

    b      = ds/(2*pi);                           % r = b*phi
    phimax = aRad/b;                              % reach r = a
    phi0   = pi;                                  % start ~half turn (centre covered)
    M      = 6000;
    phi    = linspace(phi0, phimax, M);
    dphi   = phi(2) - phi(1);

    % ---- spiral in on-sky coords (x = cross-elevation, y = elevation) ----
    x = b*phi.*cos(phi);
    y = b*phi.*sin(phi);

    % ---- map to gimbal axes: El = el0 + y ; Az = az0 + x/cos(el0) --------
    c  = 1/cos(el0Rad);                           % azimuth stretch (keyhole)
    az = c*x;                                     % offset from az0 [rad]
    el = y;                                        % offset from el0 [rad]

    azp  = gradient(az,  dphi);  elp  = gradient(el,  dphi);
    azpp = gradient(azp, dphi);  elpp = gradient(elp, dphi);
    xp   = gradient(x,   dphi);  yp   = gradient(y,   dphi);
    pskyp = sqrt(xp.^2 + yp.^2);                  % on-sky path speed / phidot

    % ---- maximum velocity curve (MVC): cap on phidot at each point -------
    eps0 = 1e-12;
    v_az  = Vaz  ./ max(abs(azp),  eps0);
    v_el  = Vel  ./ max(abs(elp),  eps0);
    v_sky = Vsky ./ max(abs(pskyp),eps0);
    a_az  = sqrt(Aaz ./ max(abs(azpp), eps0));    % pure-centripetal cap
    a_el  = sqrt(Ael ./ max(abs(elpp), eps0));
    comps = [v_az; v_el; v_sky; a_az; a_el];
    [MVC, whichLim] = min(comps, [], 1);

    % ---- velocity/PRF-only time (accel ignored) -------------------------
    vlim = min([v_az; v_el; v_sky], [], 1);
    T_vel = sum(dphi ./ max((vlim(1:end-1)+vlim(2:end))/2, 1e-12));

    % ---- forward-backward acceleration-limited pass ---------------------
    phidot = MVC;
    phidot(1) = 0;                                % start from rest at centre
    for i = 1:M-1                                 % forward: accelerate
        fd = phidot(i);
        hi = min( phiddotHi(azp(i),azpp(i),fd,Aaz), phiddotHi(elp(i),elpp(i),fd,Ael) );
        nd2 = fd^2 + 2*hi*dphi;
        phidot(i+1) = min(MVC(i+1), sqrt(max(nd2,0)));
    end
    for i = M-1:-1:1                              % backward: ensure we can brake
        fd = phidot(i+1);
        loMag = min( phiddotHi(azp(i+1),azpp(i+1),fd,Aaz), phiddotHi(elp(i+1),elpp(i+1),fd,Ael) );
        nd2 = fd^2 + 2*loMag*dphi;                % symmetric decel capability
        phidot(i) = min(phidot(i), sqrt(max(nd2,0)));
    end
    phidot = max(phidot, 1e-9);
    T_full = sum(dphi ./ ((phidot(1:end-1)+phidot(2:end))/2));

    % ---- idealised point-based metrics ----------------------------------
    K        = pi*aRad^2/ds^2 + 1;
    T_dwellLB= K / hw.PRF;
    Lsky     = sum(sqrt(diff(x).^2 + diff(y).^2));

    % ---- which constraint binds (fraction of path) ----------------------
    names = {'v_az','v_el','v_sky_PRF','a_az','a_el'};
    bf.v_az      = mean(whichLim==1);
    bf.v_el      = mean(whichLim==2);
    bf.v_sky_PRF = mean(whichLim==3);
    bf.a_az      = mean(whichLim==4);
    bf.a_el      = mean(whichLim==5);

    out.T_full   = T_full;
    out.T_vel    = T_vel;
    out.T_dwellLB= T_dwellLB;
    out.K        = K;
    out.Nturn    = aRad/ds;
    out.LskyDeg  = Lsky/deg;
    out.dsArcsec = ds/as;
    out.VskyDeg  = Vsky/deg;
    out.bindFrac = bf;
    out.bindNames= names;
    out.whichLim = whichLim;
    out.phi      = phi;
    out.phidot   = phidot;
    out.MVC      = MVC;
    out.az       = az;      % offsets [rad]
    out.el       = el;
    out.x        = x;       % on-sky offsets [rad]
    out.y        = y;
end

% -------------------------------------------------------------------------
function hi = phiddotHi(ap, app, fd, A)
% Upper bound on phiddot from |ap*phiddot + app*fd^2| <= A for one axis.
% (Symmetric: the same magnitude bounds deceleration.)
    if abs(ap) < 1e-12
        hi = 1e9;                       % axis imposes no tangential bound here
    else
        hi = (A - app*fd^2) / abs(ap);  % use |ap|: accelerate along the path
        if hi < 0, hi = 0; end
    end
end
