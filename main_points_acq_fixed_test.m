%==========================================================================
% main_points_acq_fixed_test.m
%
%   Control experiment for main_points_acq_test.m.
%
%   Question: is the acquisition failure caused by the FOU MOVING, or just by
%   its size / the scan time? To isolate this, freeze the FOU rigidly onto the
%   mean trajectory: every candidate keeps a CONSTANT angular offset from the
%   mean pointing (no 1/range breathing, no rotation, no differential motion).
%   The cloud still rides across the sky with the mean, but as a rigid pattern.
%
%   Because the frozen size depends on WHICH epoch we freeze at, two references
%   are run and compared against the moving case of Stage 4:
%     * fixed@culmination : frozen at the largest (culmination) projection
%     * fixed@start       : frozen at the small scan-start (detElev) projection
%
%   Same overlapping Archimedean spiral, same gimbal/PRF limits, same points.
%
%   Depends on the Stage-1..4 pipeline (spiralAcqSim, eciToAzElVec, ...).
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end
rng_seed = 12345;
try, rng(rng_seed); catch, rand('seed',rng_seed); randn('seed',rng_seed); end

%% ----------------------------- CONFIG ----------------------------------
station.latDeg = 36.3711; station.lonDeg = 127.3617; station.altKm = 0.10;
station.rEcef  = geodetic2ecef(station.latDeg, station.lonDeg, station.altKm);
sigmaRIC.R = 0.205; sigmaRIC.I = 0.808; sigmaRIC.C = 0.168;
passOpt.searchHours = 12; passOpt.coarseDtSec = 15; passOpt.fineDtSec = 1; passOpt.maskDeg = 10;
projOpt.nSigma = 3; projOpt.sigmaIsoRad = 0; projOpt.nEllipse = 20;

detElevDeg = 20;  Npts = 3000;  alpha = 1.0;
hw.Vaz = 20; hw.Vel = 10; hw.Aaz = 5; hw.Ael = 2; hw.PRF = 2000;
beamList = [5 10 20 50 100 200];
FOCUS = 1;
OUTDIR = 'figures'; if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

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
nSat = size(sats,1); nB = numel(beamList);
as = pi/180/3600;

L  = chol(diag([sigmaRIC.R^2, sigmaRIC.I^2, sigmaRIC.C^2]), 'lower');
dR = (L * randn(3, Npts))';

acqMove = nan(nSat,nB); acqCulm = nan(nSat,nB); acqStart = nan(nSat,nB);
Rall = cell(nSat,1);

fprintf('\n=== Fixed-FOU control | N=%d, detElev=%g, alpha=%.2f ===\n', Npts, detElevDeg, alpha);
for s = 1:nSat
    cfgBase.mode='tle'; cfgBase.tle = parseTLE(sats{s,2}, sats{s,3});
    pass = findBestPass(cfgBase, station, passOpt);
    cfg = cfgBase; cfg.tSec = pass.tSec; trk = generateTrack(cfg);
    cov = generateCovariance(trk, sigmaRIC);
    fou = projectFoU(cov.SigmaEci, trk.rEci, trk.jd, station, projOpt);

    win = find(fou.altDeg >= detElevDeg); i0 = win(1); i1 = win(end);
    idx = i0:i1; M = numel(idx);
    tGrid = trk.tSec(idx); meanEl = deg2rad(fou.altDeg(idx));
    meanAz = deg2rad(fou.azDeg(idx));
    [~, ipk] = max(meanEl);

    dCross = zeros(Npts,M); dEl = zeros(Npts,M);
    for jj = 1:M
        m = idx(jj);
        Q = [cov.Rhat(:,m), cov.Ihat(:,m), cov.Chat(:,m)];
        candEci = trk.rEci(:,m) + Q*dR.';
        [caz, cel] = eciToAzElVec(candEci, gmstRad(trk.jd(m)), station);
        dA = mod(caz - meanAz(jj) + pi, 2*pi) - pi;
        dCross(:,jj) = (dA .* cos(meanEl(jj))).';
        dEl(:,jj)    = (cel - meanEl(jj)).';
    end
    rho = sqrt(dCross.^2 + dEl.^2);

    % frozen (rigid co-moving) versions: constant local offset over time
    dCrossC = repmat(dCross(:,ipk), 1, M); dElC = repmat(dEl(:,ipk), 1, M);
    rhoC = sqrt(dCrossC.^2 + dElC.^2);
    dCrossS = repmat(dCross(:,1),   1, M); dElS = repmat(dEl(:,1),   1, M);
    rhoS = sqrt(dCrossS.^2 + dElS.^2);

    RmaxMove = 1.15*max(rho(:));  RmaxC = 1.15*max(rhoC(:));  RmaxS = 1.15*max(rhoS(:));

    fprintf('\n%-16s  window %.1f min, a: %.0f(start)->%.0f(culm) urad\n', ...
        sats{s,1}, (tGrid(end)-tGrid(1))/60, fou.aRad(i0)*1e6, fou.aRad(idx(ipk))*1e6);
    fprintf('   %6s %9s %13s %11s\n','beam"','moving','fixed@culm','fixed@start');
    for bi = 1:nB
        Rbeam = 0.5*beamList(bi)*as;
        aM = spiralAcqSim(tGrid, meanEl, dCross,  dEl,  rho,  Rbeam, alpha, hw, RmaxMove).acqRate;
        aC = spiralAcqSim(tGrid, meanEl, dCrossC, dElC, rhoC, Rbeam, alpha, hw, RmaxC).acqRate;
        aS = spiralAcqSim(tGrid, meanEl, dCrossS, dElS, rhoS, Rbeam, alpha, hw, RmaxS).acqRate;
        acqMove(s,bi)=aM; acqCulm(s,bi)=aC; acqStart(s,bi)=aS;
        fprintf('   %6d %8.1f%% %12.1f%% %10.1f%%\n', beamList(bi), 100*aM, 100*aC, 100*aS);
    end
    Rall{s} = struct('name',sats{s,1});
end

%% --------------------------- PLOT --------------------------------------
s = FOCUS;
fig = figure('Position',[80 80 820 560],'Color','w','Visible','off');
hold on; grid on; box on;
plot(beamList, 100*acqMove(s,:),  '-o', 'Color',[0.85 0.15 0.10], 'LineWidth',1.9, 'MarkerFaceColor',[0.85 0.15 0.10]);
plot(beamList, 100*acqCulm(s,:),  '-s', 'Color',[0.20 0.45 0.85], 'LineWidth',1.9, 'MarkerFaceColor',[0.20 0.45 0.85]);
plot(beamList, 100*acqStart(s,:), '-^', 'Color',[0.30 0.65 0.30], 'LineWidth',1.9, 'MarkerFaceColor',[0.30 0.65 0.30]);
set(gca,'XScale','log'); ylim([0 105]);
xlabel('beam width [arcsec]'); ylabel('acquisition success rate [%]');
title(sprintf('%s: moving FOU vs FOU frozen onto the trajectory', Rall{s}.name));
legend({'moving FOU (Stage 4)','fixed @ culmination size','fixed @ start size'}, ...
       'Location','southeast'); legend boxoff;
print(fig, fullfile(OUTDIR,'acq_fixed_vs_moving.png'), '-dpng','-r120'); close(fig);

fprintf('\nFigure: figures/acq_fixed_vs_moving.png\n');
fprintf('Done.\n');
