%==========================================================================
% main_sept_spiral_test.m
%
%   Particle-based conceptual study of a CONVENTIONAL Archimedean-spiral scan,
%   sweeping the beam control frequency f in {100,250,500,1000,2000} Hz.
%
%   Spiral: overlap = 1  (alpha = 1  =>  R_s = 2*Rbeam, adjacent footprints
%   exactly tangent, so interstitial gaps remain -> even a FIXED FOU is not
%   fully harvested). It is drawn for the LARGEST FOU (satellite near zenith,
%   TLE-based) out to the FOU major axis a. Scan point k is visited at
%   t_k = t_start + (k-1)/f  (pure control-frequency timing; no gimbal
%   dynamics). Scan ends when the spiral passes the major axis (K points).
%
%   Metric: harvest probability = fraction of N Gaussian particles (candidate
%   satellite positions) covered within Rbeam of some scan point at that
%   point's visit time.
%
%   (A) STATIONARY : FOU and particles frozen at the zenith epoch. Harvest is
%       frequency-independent (same scan points); only the scan TIME = K/f
%       changes (inversely with f).
%   (B) ORBITING   : same spiral / same scan time K/f, spiral centred on the
%       moving mean orbit, particles moving. Two starts:
%         (1) at first detectable elevation (10 deg),
%         (2) at the largest FOU (culmination ~ zenith).
%       Harvest is computed at scan end and compared with (A).
%
%   Files: buildSpiralPoints.m, harvestScan.m, precomputeOffsets.m + pipeline.
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end
rng_seed = 12345; try, rng(rng_seed); catch, randn('seed',rng_seed); end

%% ----------------------------- CONFIG ----------------------------------
station.latDeg = 36.3711; station.lonDeg = 127.3617; station.altKm = 0.10;
station.rEcef  = geodetic2ecef(station.latDeg, station.lonDeg, station.altKm);
sigmaRIC.R = 0.205; sigmaRIC.I = 0.808; sigmaRIC.C = 0.168;

passOpt.searchHours = 12; passOpt.coarseDtSec = 15; passOpt.fineDtSec = 2; passOpt.maskDeg = 10;
projOpt.nSigma = 3; projOpt.sigmaIsoRad = 0; projOpt.nEllipse = 12;

Npts       = 40000;
beamArcsec = 50;              % beam width (footprint diameter) [arcsec]
alpha      = 1.0;            % overlap = 1  (tangent footprints, gaps remain)
fList      = [100 250 500 1000 2000];       % beam control frequency [Hz]
detElevDeg = 10;
dtFine     = 0.1;            % fine grid for moving-particle offsets [s]
OUTDIR = 'figures'; if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

tle = parseTLE('1 66601U 25270P   26238.87143828  .00021431  00000-0  78601-3 0  9993', ...
               '2 66601  43.0018 175.0289 0001887 251.0337 109.0311 15.27581436 43570');
as = pi/180/3600; Rbeam = 0.5*beamArcsec*as; nF = numel(fList);

%% ------------------- pass geometry: peak (zenith) FOU ------------------
cfgBase.mode='tle'; cfgBase.tle=tle;
pass = findBestPass(cfgBase, station, passOpt);
cfg = cfgBase; cfg.tSec = pass.tSec; trk = generateTrack(cfg);
cov = generateCovariance(trk, sigmaRIC);
fou = projectFoU(cov.SigmaEci, trk.rEci, trk.jd, station, projOpt);
[~, ipk] = max(fou.altDeg);
i10 = find(fou.altDeg >= detElevDeg, 1, 'first');
aZen  = fou.aRad(ipk);                 % major axis of the largest (zenith) FOU
tPeak = trk.tSec(ipk);  t10 = trk.tSec(i10);

%% ------------------- spiral (overlap=1) to the major axis --------------
[sx, sy, K, Rs] = buildSpiralPoints(aZen, Rbeam, alpha);
tScan = K ./ fList;                     % scan-completion time per frequency [s]

fprintf('\n=== Sept spiral study | Starlink 35978 | beam=%g", alpha=%g (overlap=1) ===\n', beamArcsec, alpha);
fprintf('peak elev %.1f deg, a_zenith=%.0f urad, R_s=%.0f urad, K=%d scan points\n', ...
        fou.altDeg(ipk), aZen*1e6, Rs*1e6, K);

%% ------------------- particles (Gaussian) -----------------------------
L  = chol(diag([sigmaRIC.R^2, sigmaRIC.I^2, sigmaRIC.C^2]), 'lower');
dR = (L * randn(3, Npts)).';            % Npts x 3

%% ================== (A) STATIONARY =====================================
[sX, sY] = precomputeOffsets(tPeak, tle, station, sigmaRIC, dR);   % frozen at peak
harvStat = harvestScan(sx, sy, zeros(1,K), sX, sY, 0, Rbeam);
Hstat = mean(harvStat);
fprintf('\n(A) STATIONARY (frozen zenith FOU):  harvest = %.1f%%  (freq-independent)\n', 100*Hstat);
fprintf('    scan time = K/f:');
for i=1:nF, fprintf('  %gHz->%.2fs', fList(i), tScan(i)); end
fprintf('\n');

%% ================== (B) ORBITING =======================================
Tmax = max(tScan);                       % longest scan (lowest frequency)
t1 = t10  :dtFine: t10 +Tmax+dtFine;     % window for start@10deg
t2 = tPeak:dtFine: tPeak+Tmax+dtFine;    % window for start@zenith
[o1X,o1Y] = precomputeOffsets(t1, tle, station, sigmaRIC, dR);
[o2X,o2Y] = precomputeOffsets(t2, tle, station, sigmaRIC, dR);

H1 = zeros(1,nF); H2 = zeros(1,nF);
fprintf('\n(B) ORBITING (same scan time K/f, spiral centred on mean orbit):\n');
fprintf('   %8s | %-12s %-12s | %-12s\n','f[Hz]','harvest@10deg','harvest@zenith','scan time');
for i = 1:nF
    tk = t10  + (0:K-1)/fList(i);   H1(i) = mean(harvestScan(sx,sy,tk, o1X,o1Y,t1, Rbeam));
    tk = tPeak+ (0:K-1)/fList(i);   H2(i) = mean(harvestScan(sx,sy,tk, o2X,o2Y,t2, Rbeam));
    fprintf('   %8d | %10.1f%% %12.1f%% | %8.2f s\n', fList(i), 100*H1(i), 100*H2(i), tScan(i));
end

%% ----------------------------- PLOTS -----------------------------------
% (1) harvest vs frequency
fig = figure('Position',[70 70 760 560],'Color','w','Visible','off');
hold on; grid on; box on;
plot(fList, 100*Hstat*ones(1,nF), 'k--', 'LineWidth',1.6);
plot(fList, 100*H1, '-o', 'Color',[0.20 0.55 0.30], 'LineWidth',1.9,'MarkerFaceColor',[0.20 0.55 0.30]);
plot(fList, 100*H2, '-s', 'Color',[0.85 0.20 0.15], 'LineWidth',1.9,'MarkerFaceColor',[0.85 0.20 0.15]);
set(gca,'XScale','log'); ylim([0 103]);
xlabel('beam control frequency [Hz]'); ylabel('harvest probability [%]');
title(sprintf('Conventional spiral (overlap=1, beam %g arcsec): harvest vs frequency', beamArcsec));
legend({'stationary (frozen FOU)','orbiting: start 10 deg','orbiting: start zenith'}, ...
       'Location','southeast'); legend boxoff;
print(fig, fullfile(OUTDIR,'sept_spiral_harvest.png'), '-dpng','-r120'); close(fig);

% (2) scan time vs frequency
fig = figure('Position',[70 70 620 460],'Color','w','Visible','off');
loglog(fList, tScan, '-o', 'Color',[0.2 0.2 0.2], 'LineWidth',1.8,'MarkerFaceColor',[0.2 0.2 0.2]);
grid on; box on; xlabel('beam control frequency [Hz]'); ylabel('scan-completion time K/f [s]');
title(sprintf('Scan time (K=%d points): time = K/f', K));
print(fig, fullfile(OUTDIR,'sept_spiral_time.png'), '-dpng','-r120'); close(fig);

% (3) scan-point placement + beam footprints over the STATIONARY FOU
urad = 1e6; th = linspace(0,2*pi,72);
[V,D] = eig(fou.SigmaAng(:,:,ipk)); lam = max(diag(D),0);
ex = projOpt.nSigma*(V(1,1)*sqrt(lam(1))*cos(th) + V(1,2)*sqrt(lam(2))*sin(th));
ey = projOpt.nSigma*(V(2,1)*sqrt(lam(1))*cos(th) + V(2,2)*sqrt(lam(2))*sin(th));
Rb_u = Rbeam*urad; cc = cos(th); ss = sin(th);
ac = harvStat; mi = ~harvStat;

fig = figure('Position',[50 50 1240 580],'Color','w','Visible','off');

% full view: FOU ellipse + all scan points + harvested/missed particles
subplot(1,2,1); hold on; grid on; box on; axis equal;
plot(sX(ac)*urad, sY(ac)*urad, '.', 'Color',[0.25 0.60 0.30], 'MarkerSize',2);
plot(sX(mi)*urad, sY(mi)*urad, '.', 'Color',[0.85 0.15 0.10], 'MarkerSize',3);
plot(sx*urad, sy*urad, '.', 'Color',[0.15 0.30 0.85], 'MarkerSize',3);
plot(ex*urad, ey*urad, 'k--', 'LineWidth',2);
xlabel('cross-elevation [urad]'); ylabel('elevation [urad]');
title(sprintf('spiral over 3-sigma FOU  (harvest %.1f%%, K=%d)', 100*Hstat, K));
legend({'harvested','missed','scan points','FOU 3-sigma'}, ...
       'Location','southoutside','Orientation','horizontal'); legend boxoff;

% zoom: beam footprints (overlap=1) showing tangency + interstitial gaps
subplot(1,2,2); hold on; grid on; box on; axis equal;
zl = 1200;                                        % zoom half-width [urad]
in = find(sqrt(sx.^2+sy.^2)*urad <= zl+Rb_u);
for q = in
    fill(sx(q)*urad+Rb_u*cc, sy(q)*urad+Rb_u*ss, [0.55 0.70 1.0], ...
         'FaceAlpha',0.25, 'EdgeColor',[0.45 0.55 0.9], 'LineWidth',0.3);
end
plot(sx(in)*urad, sy(in)*urad, 'k.', 'MarkerSize',4);
pin = (abs(sX*urad)<=zl) & (abs(sY*urad)<=zl);
plot(sX(pin & mi)*urad, sY(pin & mi)*urad, '.', 'Color',[0.85 0.10 0.05], 'MarkerSize',6);
plot(sX(pin & ac)*urad, sY(pin & ac)*urad, '.', 'Color',[0.10 0.45 0.15], 'MarkerSize',3);
xlim([-zl zl]); ylim([-zl zl]);
xlabel('cross-elevation [urad]'); ylabel('elevation [urad]');
title(sprintf('beam footprints (overlap=1, Rbeam=%.0f urad): tangent, gaps -> misses', Rb_u));
print(fig, fullfile(OUTDIR,'sept_spiral_layout.png'), '-dpng','-r120'); close(fig);

% (4) ORBITING: scanned vs unscanned particles, plotted in the SAME zenith
%     projection frame as the stationary case. The effective spiral changes as
%     the FOU rescales/rotates with elevation, so footprint overlap is not
%     meaningful here -> only the scanned/unscanned outcome is shown. Each
%     particle sits at its zenith-frame offset (sX,sY); colour = orbiting result.
fRep = fList(1);                                             % representative freq (longest scan, most motion)
hO1 = harvestScan(sx,sy, t10  +(0:K-1)/fRep, o1X,o1Y,t1, Rbeam);   % start @ 10 deg
hO2 = harvestScan(sx,sy, tPeak+(0:K-1)/fRep, o2X,o2Y,t2, Rbeam);   % start @ zenith

fig = figure('Position',[50 50 1240 600],'Color','w','Visible','off');
subplot(1,2,1); hold on; grid on; box on; axis equal;
plot(sX(hO1) *urad, sY(hO1) *urad, '.', 'Color',[0.25 0.60 0.30], 'MarkerSize',2);
plot(sX(~hO1)*urad, sY(~hO1)*urad, '.', 'Color',[0.85 0.15 0.10], 'MarkerSize',3);
plot(ex*urad, ey*urad, 'k--', 'LineWidth',1.6);
xlabel('cross-elevation [urad]'); ylabel('elevation [urad]');
title(sprintf('orbiting, start @ 10 deg:  scanned %.1f%%  (stationary %.1f%%)', 100*mean(hO1), 100*Hstat));
legend({'scanned','unscanned','zenith 3-sigma FOU'}, ...
       'Location','southoutside','Orientation','horizontal'); legend boxoff;

subplot(1,2,2); hold on; grid on; box on; axis equal;
plot(sX(hO2) *urad, sY(hO2) *urad, '.', 'Color',[0.25 0.60 0.30], 'MarkerSize',2);
plot(sX(~hO2)*urad, sY(~hO2)*urad, '.', 'Color',[0.85 0.15 0.10], 'MarkerSize',3);
plot(ex*urad, ey*urad, 'k--', 'LineWidth',1.6);
xlabel('cross-elevation [urad]'); ylabel('elevation [urad]');
title(sprintf('orbiting, start @ zenith:  scanned %.1f%%  (stationary %.1f%%)', 100*mean(hO2), 100*Hstat));
legend({'scanned','unscanned','zenith 3-sigma FOU'}, ...
       'Location','southoutside','Orientation','horizontal'); legend boxoff;
print(fig, fullfile(OUTDIR,'sept_spiral_orbit_layout.png'), '-dpng','-r120'); close(fig);

fprintf('\nFigures: figures/sept_spiral_harvest.png, figures/sept_spiral_time.png,\n');
fprintf('         figures/sept_spiral_layout.png, figures/sept_spiral_orbit_layout.png\n');
fprintf('Done.\n');
