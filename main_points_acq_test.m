%==========================================================================
% main_points_acq_test.m
%
%   Monte-Carlo acquisition test with MOVING candidate points.
%
%   Instead of an analytic FOU, draw N (>=3000) candidate points from the
%   in-track-elongated orbit covariance. Each candidate is a fixed RIC offset
%   (a hypothesis about the true orbit); propagated through ECI->ENU it traces
%   its own moving path on the sky, so the whole cloud translates, rotates and
%   breathes (1/range) along the pass -- exactly the analytic FOU, but sampled.
%
%   Above a detectable elevation the gimbal scans the mean path with an
%   OVERLAPPING Archimedean spiral (points at theta_k = sqrt(4*pi*k),
%   r_k = p*sqrt(k/pi), pitch p = alpha*Rbeam, alpha<2 => overlap). Each
%   candidate is "acquired" the first time a scan point lands within Rbeam of
%   the candidate's CURRENT position at that point's visit time.
%
%   Purpose: show that because the FOU MOVES while the spiral scans, some
%   candidates are never coincident with the beam and the scan can FAIL
%   (acquisition rate < 100%), especially for fine beams / high elevations.
%
%   Depends on the Stage-1..3 pipeline + spiralAcqSim.m + eciToAzElVec.m.
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end
rng_seed = 12345;
try, rng(rng_seed); catch, rand('seed',rng_seed); randn('seed',rng_seed); end

%% ----------------------------- CONFIG ----------------------------------
station.latDeg = 36.3711; station.lonDeg = 127.3617; station.altKm = 0.10;
station.rEcef  = geodetic2ecef(station.latDeg, station.lonDeg, station.altKm);

sigmaRIC.R = 0.205; sigmaRIC.I = 0.808; sigmaRIC.C = 0.168;   % RIC 1-sigma [km]

passOpt.searchHours = 12; passOpt.coarseDtSec = 15; passOpt.fineDtSec = 1; passOpt.maskDeg = 10;
projOpt.nSigma = 3; projOpt.sigmaIsoRad = 0; projOpt.nEllipse = 20;

detElevDeg = 20;              % detectable elevation: scan only above this
Npts       = 3000;           % number of covariance sample points
alpha      = 1.0;            % spiral overlap/spacing factor (pitch = alpha*Rbeam)

hw.Vaz = 20; hw.Vel = 10; hw.Aaz = 5; hw.Ael = 2; hw.PRF = 2000;

beamList = [5 10 20 50 100 200];    % beam width [arcsec]
FOCUS = 1; FOCUS_BEAM = 10;         % satellite / beam for detail plots
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
nSat = size(sats,1);
as = pi/180/3600;

%% --------------------------- RUN ---------------------------------------
fprintf('\n=== Points acquisition test | N=%d, detElev=%g deg, alpha=%.2f, PRF=%d ===\n', ...
        Npts, detElevDeg, alpha, hw.PRF);

% covariance sample offsets in RIC [km], shared across satellites
L  = chol(diag([sigmaRIC.R^2, sigmaRIC.I^2, sigmaRIC.C^2]), 'lower');
dR = (L * randn(3, Npts))';               % Npts x 3

Rall = cell(nSat,1);
for s = 1:nSat
    cfgBase.mode='tle'; cfgBase.tle = parseTLE(sats{s,2}, sats{s,3});
    pass = findBestPass(cfgBase, station, passOpt);
    cfg = cfgBase; cfg.tSec = pass.tSec; trk = generateTrack(cfg);
    cov = generateCovariance(trk, sigmaRIC);
    fou = projectFoU(cov.SigmaEci, trk.rEci, trk.jd, station, projOpt);

    % detectable window (contiguous block above detElevDeg)
    win = find(fou.altDeg >= detElevDeg);
    i0 = win(1); i1 = win(end);
    idx = i0:i1;  M = numel(idx);
    tGrid  = trk.tSec(idx);
    meanAz = deg2rad(fou.azDeg(idx));
    meanEl = deg2rad(fou.altDeg(idx));
    [~, ipkW] = max(meanEl);

    % candidate local offsets on the window grid
    dCross = zeros(Npts, M); dEl = zeros(Npts, M); rho = zeros(Npts, M);
    for jj = 1:M
        m = idx(jj);
        Q = [cov.Rhat(:,m), cov.Ihat(:,m), cov.Chat(:,m)];   % RIC->ECI
        candEci = trk.rEci(:,m) + Q * dR.';                  % 3 x N
        [caz, cel] = eciToAzElVec(candEci, gmstRad(trk.jd(m)), station);
        dA = mod(caz - meanAz(jj) + pi, 2*pi) - pi;          % wrapped az diff
        dCross(:,jj) = (dA .* cos(meanEl(jj))).';
        dEl(:,jj)    = (cel - meanEl(jj)).';
    end
    rho = sqrt(dCross.^2 + dEl.^2);
    RmaxCover = 1.15 * max(rho(:));

    S = struct('name',sats{s,1},'tGrid',tGrid,'meanEl',meanEl,'ipk',ipkW, ...
               'dCross',dCross,'dEl',dEl,'rho',rho,'a',fou.aRad(idx));
    S.acqRate = zeros(1,numel(beamList)); S.Kpts = zeros(1,numel(beamList));
    S.reached = false(1,numel(beamList)); S.sims = cell(1,numel(beamList));

    fprintf('\n%-16s  scan window %.1f min (elev %.0f->%.0f->%.0f deg), FOU a max %.0f urad\n', ...
        sats{s,1}, (tGrid(end)-tGrid(1))/60, rad2deg(meanEl(1)), rad2deg(meanEl(ipkW)), ...
        rad2deg(meanEl(end)), max(fou.aRad(idx))*1e6);
    fprintf('   %6s %8s %10s %9s  %s\n', 'beam"','ds"','acq rate','K pts','spiral reached FOU edge?');
    for bi = 1:numel(beamList)
        Rbeam = 0.5*beamList(bi)*as;
        sim = spiralAcqSim(tGrid, meanEl, dCross, dEl, rho, Rbeam, alpha, hw, RmaxCover);
        S.acqRate(bi) = sim.acqRate; S.Kpts(bi) = sim.Kpts; S.reached(bi) = sim.reached;
        S.sims{bi} = sim;
        fprintf('   %6d %8.1f %9.1f%% %9d  %s\n', beamList(bi), alpha*beamList(bi), ...
            100*sim.acqRate, sim.Kpts, ternary(sim.reached,'yes','NO (ran out of pass)'));
    end
    Rall{s} = S;
end

%% --------------------------- PLOTS -------------------------------------
% (1) acquisition rate vs beam width
fig = figure('Position',[80 80 760 560],'Color','w','Visible','off');
hold on; grid on; box on; cols = lines(nSat);
for s = 1:nSat
    plot(beamList, 100*Rall{s}.acqRate, '-o', 'Color',cols(s,:), 'LineWidth',1.8, 'MarkerFaceColor',cols(s,:));
end
set(gca,'XScale','log'); ylim([0 105]);
xlabel('beam width [arcsec]'); ylabel('acquisition success rate [%]');
title(sprintf('Acquisition of %d moving points vs beam width (detElev %g deg)', Npts, detElevDeg));
legend(cellfun(@(r) r.name, Rall,'UniformOutput',false),'Location','southeast'); legend boxoff;
print(fig, fullfile(OUTDIR,'acq_rate_vs_beam.png'), '-dpng','-r120'); close(fig);

% (2) sky-plane snapshot at culmination: acquired vs missed + scan points
S = Rall{FOCUS}; bi = find(beamList==FOCUS_BEAM,1); sim = S.sims{bi}; urad = 1e6; ipk = S.ipk;
fig = figure('Position',[80 80 720 680],'Color','w','Visible','off');
hold on; grid on; box on; axis equal;
plot(sim.sx*urad, sim.sy*urad, '.', 'Color',[0.7 0.7 0.7], 'MarkerSize',3);
mi = ~sim.acquired; ac = sim.acquired;
plot(S.dCross(ac,ipk)*urad, S.dEl(ac,ipk)*urad, '.', 'Color',[0.20 0.55 0.25], 'MarkerSize',5);
plot(S.dCross(mi,ipk)*urad, S.dEl(mi,ipk)*urad, '.', 'Color',[0.85 0.15 0.10], 'MarkerSize',6);
xlabel('cross-elevation [urad]'); ylabel('elevation [urad]');
title(sprintf('%s, %d arcsec: candidates at culmination (acq %.1f%%)', S.name, FOCUS_BEAM, 100*sim.acqRate));
legend({'scan points','acquired','missed'},'Location','southoutside','Orientation','horizontal'); legend boxoff;
print(fig, fullfile(OUTDIR,'acq_sky_snapshot.png'), '-dpng','-r120'); close(fig);

% (3) cumulative acquisition vs time, focus satellite, several beams
fig = figure('Position',[80 80 780 560],'Color','w','Visible','off');
hold on; grid on; box on; cmap = jet(numel(beamList)); labs = {};
trel = S.tGrid - S.tGrid(1); tcul = trel(ipk);
for bi = 1:numel(beamList)
    th = sort(S.sims{bi}.tHit(~isnan(S.sims{bi}.tHit)));
    if isempty(th), continue; end
    frac = 100*(1:numel(th))/Npts;
    stairs([0; th(:)], [0; frac(:)], '-', 'Color',cmap(bi,:), 'LineWidth',1.5);
    labs{end+1} = sprintf('%d arcsec (%.0f%%)', beamList(bi), 100*S.acqRate(bi)); %#ok<AGROW>
end
yl = ylim; plot(tcul*[1 1], [0 105], ':', 'Color',[0.4 0.4 0.4]);
text(tcul, 5, ' culmination','Color',[0.4 0.4 0.4]);
xlabel('time since scan start [s]'); ylabel('cumulative acquired [%]'); ylim([0 105]);
title(sprintf('%s: when candidates are acquired', S.name));
legend(labs,'Location','southeast'); legend boxoff;
print(fig, fullfile(OUTDIR,'acq_cumulative.png'), '-dpng','-r120'); close(fig);

fprintf('\nFigures: acq_rate_vs_beam.png, acq_sky_snapshot.png, acq_cumulative.png\n');
fprintf('Done.\n');
