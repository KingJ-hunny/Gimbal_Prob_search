%==========================================================================
% main_projection_test.m
%
%   Demonstration / sanity test for the FOU (Field of Uncertainty) pipeline.
%
%   Goal: show how the along-track-elongated orbit-prediction covariance,
%   propagated along a satellite pass and projected through the ECI->ECEF->ENU
%   chain onto the plane of sky, produces an ELLIPTICAL angular FOU whose
%   size and orientation change with orbit type, geometry and motion:
%
%     * angular size  ~ sigma_pos / range   -> shrinks at low elevation
%       (far), grows near culmination (near) — the 1/range law of Kim et al.
%     * orientation   follows the in-track direction painted on the sky,
%       so the ellipse rotates as the ground track curves overhead.
%
%   Ground station is fixed in ECEF (Korea). Realistic orbits come from TLEs
%   propagated with SGP4 (sgp4.m). An analytic Kepler mode is also wired in
%   (set USE_TLE=false) so the same test can be run for MEO/GEO/HEO later.
%
%   Files used:
%     generateTrack.m  (Fn 1)  parseTLE.m  sgp4.m
%     generateCovariance.m (Fn 2)
%     projectFoU.m  findBestPass.m  computeAzAlt.m + frame helpers
%==========================================================================
clear; clc; close all;

% Octave: use the gnuplot toolkit for headless PNG export (no-op in MATLAB).
if exist('OCTAVE_VERSION','builtin')
    try, graphics_toolkit('gnuplot'); catch, end
end

%% ----------------------------- CONFIG ----------------------------------
% Ground station (Korea) — Daejeon / KAIST. Edit freely.
station.latDeg = 36.3711;
station.lonDeg = 127.3617;
station.altKm  = 0.10;
station.rEcef  = geodetic2ecef(station.latDeg, station.lonDeg, station.altKm);

% Orbit-prediction error, 1-sigma in RIC [km] (in-track dominated).
% Defaults ~ Kim et al. (Optics Comm. 620, 2026) Table 2, 24-48 h TLE age.
sigmaRIC.R = 0.205;    % radial
sigmaRIC.I = 0.808;    % in-track   (dominant -> ellipse major axis)
sigmaRIC.C = 0.168;    % cross-track

% Projection options
projOpt.nSigma      = 3;          % FOU drawn at 3-sigma (Kim et al.)
projOpt.sigmaIsoRad = 0;          % isotropic angular term (attitude/boresight); 0 = pure orbit
projOpt.nEllipse    = 120;

% Pass search options
passOpt.searchHours = 12;
passOpt.coarseDtSec = 15;
passOpt.fineDtSec   = 2;
passOpt.maskDeg     = 10;         % horizon mask

DO_PLOTS = true;
OUTDIR   = 'figures';

% -------- Realistic (TLE) targets --------------------------------------
USE_TLE = true;
sats = {
  'Starlink 35978', ...
  '1 66601U 25270P   26238.87143828  .00021431  00000-0  78601-3 0  9993', ...
  '2 66601  43.0018 175.0289 0001887 251.0337 109.0311 15.27581436 43570';
  'Starlink 31742', ...
  '1 59499U 24071R   26238.99246617 -.00000315  00000-0 -19221-5 0  9996', ...
  '2 59499  43.0020  97.7513 0001331 249.4654 110.6055 15.27573182134255';
  'COSMOS 1052',    ...
  '1 11129U 78109B   26238.92695670  .00000001  00000-0  10847-3 0  9998', ...
  '2 11129  74.0149  56.5202 0049285  42.3140  28.8688 12.53730015184182';
};

% -------- Analytic (Kepler) targets, used when USE_TLE=false -----------
% (kept here so the same test can exercise MEO/GEO/HEO via generateTrack)
keplerSats = {
  'MEO (a=20200)',  struct('aKm',20200,'ecc',0.001,'incDeg',55,'raanDeg',80, 'argpDeg',0,  'M0Deg',0,   'epochJD',2461279.5);
  'GEO (a=42164)',  struct('aKm',42164,'ecc',0.0002,'incDeg',0.05,'raanDeg',0,'argpDeg',0, 'M0Deg',120, 'epochJD',2461279.5);
  'Molniya HEO',    struct('aKm',26600,'ecc',0.74,'incDeg',63.4,'raanDeg',90,'argpDeg',270,'M0Deg',10,  'epochJD',2461279.5);
};

%% --------------------------- RUN LOOP ----------------------------------
if USE_TLE, nSat = size(sats,1); else, nSat = size(keplerSats,1); end
R = cell(nSat,1);

fprintf('\n=== FOU projection test | station lat %.3f lon %.3f ===\n', ...
        station.latDeg, station.lonDeg);
fprintf('RIC 1-sigma [km]: R=%.3f  I=%.3f  C=%.3f   (k_orb = %.3f)\n\n', ...
        sigmaRIC.R, sigmaRIC.I, sigmaRIC.C, min(sigmaRIC.C,sigmaRIC.R)/sigmaRIC.I);

for s = 1:nSat
    % ---- build the base config for this target ----
    cfgBase = struct();
    if USE_TLE
        name  = sats{s,1};
        cfgBase.mode = 'tle';
        cfgBase.tle  = parseTLE(sats{s,2}, sats{s,3});
    else
        name  = keplerSats{s,1};
        cfgBase.mode     = 'kepler';
        cfgBase.elements = keplerSats{s,2};
        cfgBase.useJ2    = true;
    end

    % ---- find the best pass over the station ----
    pass = findBestPass(cfgBase, station, passOpt);
    if ~pass.found
        fprintf('%-16s : no pass above %g deg within %g h (peak alt %.1f deg)\n', ...
                name, passOpt.maskDeg, passOpt.searchHours, pass.peakAltDeg);
        R{s} = struct('name',name,'found',false);
        continue;
    end

    % ---- generate the fine track over the pass ----
    cfg = cfgBase; cfg.tSec = pass.tSec;
    trk = generateTrack(cfg);

    % ---- (Fn 2) covariance, then project to angular FOU ----
    cov = generateCovariance(trk, sigmaRIC);
    fou = projectFoU(cov.SigmaEci, trk.rEci, trk.jd, station, projOpt);

    % ---- store & report ----
    [~, ipk] = max(fou.altDeg);
    urad = 1e6;
    res = struct('name',name,'found',true,'trk',trk,'cov',cov,'fou',fou, ...
                 'pass',pass,'ipk',ipk,'label',trk.label);
    R{s} = res;

    fprintf('%-16s [%s]\n', name, trk.label);
    fprintf('   pass: dur %.1f min, peak elev %.1f deg\n', pass.durationSec/60, pass.peakAltDeg);
    fprintf('   range      : %.0f km (rise) -> %.0f km (peak) -> %.0f km (set)\n', ...
            fou.rangeKm(1), fou.rangeKm(ipk), fou.rangeKm(end));
    fprintf('   FOU major a (3sig): %.1f urad (rise) -> %.1f urad (peak) -> %.1f urad (set)\n', ...
            fou.aRad(1)*urad, fou.aRad(ipk)*urad, fou.aRad(end)*urad);
    fprintf('   FOU minor b (3sig): %.1f urad (peak)\n', fou.bRad(ipk)*urad);
    fprintf('   axis ratio k=b/a  : %.3f (rise) -> %.3f (peak) -> %.3f (set)\n', ...
            fou.k(1), fou.k(ipk), fou.k(end));
    fprintf('   ellipse tilt      : %.1f deg (rise) -> %.1f deg (peak) -> %.1f deg (set)\n\n', ...
            fou.angleDeg(1), fou.angleDeg(ipk), fou.angleDeg(end));
end

%% ----------------------------- PLOTS -----------------------------------
if DO_PLOTS
    if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end
    for s = 1:nSat
        if isempty(R{s}) || ~R{s}.found, continue; end
        plotSatelliteFOU(R{s}, projOpt, OUTDIR);
    end
    plotSummary(R, OUTDIR);
    fprintf('Figures written to ./%s/\n', OUTDIR);
end

fprintf('Done.\n');
