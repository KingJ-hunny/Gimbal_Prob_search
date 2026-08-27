%==========================================================================
% main_points_acq_test.m
%
%   Monte-Carlo acquisition test with MOVING candidate points, using the
%   constant-linear-velocity Archimedean-spiral scan-point allocation:
%
%       r(theta) = theta * R_s/(2*pi),   R_s = 2*alpha*Rbeam,
%       theta_k  = sqrt(4*pi*k)      =>   r_k = R_s*sqrt(k/pi),
%
%   so one turn advances r by exactly R_s and the along-arc point spacing is
%   R_s. alpha is the overlap factor:
%       alpha = 1/sqrt(2)  -> R_s = sqrt(2)*Rbeam : full (overlapping) coverage
%       alpha = 1          -> R_s = 2*Rbeam       : footprints tangent (gaps)
%
%   N candidate points are drawn from the in-track covariance; each is a fixed
%   RIC offset propagated through ECI->ENU, so the cloud reproduces the moving
%   analytic FOU. A candidate is acquired the first time a scan point lands
%   within Rbeam of its CURRENT position at that point's visit time.
%
%   This run: N = 40000, both alpha values, all beam widths, 3 satellites.
%   Reports the overall acquisition rate AND the rim acquisition rate (the
%   outer 10% of the cloud by culmination radius) -- the rim is where the
%   Gaussian tail lives and where the overall metric hides the misses.
%
%   Depends on the Stage-1..4 pipeline + spiralAcqSim.m + eciToAzElVec.m.
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end
rng_seed = 12345;
try, rng(rng_seed); catch, rand('seed',rng_seed); randn('seed',rng_seed); end

%% ----------------------------- CONFIG ----------------------------------
station.latDeg = 36.3711; station.lonDeg = 127.3617; station.altKm = 0.10;
station.rEcef  = geodetic2ecef(station.latDeg, station.lonDeg, station.altKm);
sigmaRIC.R = 0.205; sigmaRIC.I = 0.808; sigmaRIC.C = 0.168;   % RIC 1-sigma [km]

passOpt.searchHours = 12; passOpt.coarseDtSec = 15; passOpt.fineDtSec = 2; passOpt.maskDeg = 10;
projOpt.nSigma = 3; projOpt.sigmaIsoRad = 0; projOpt.nEllipse = 12;

detElevDeg = 20;
Npts       = 40000;
alphaList  = [1/sqrt(2), 1.0];        % overlap factor (R_s = 2*alpha*Rbeam)
rimFrac    = 0.10;                     % outer fraction defining the "rim"

hw.Vaz = 20; hw.Vel = 10; hw.Aaz = 5; hw.Ael = 2; hw.PRF = 2000;
beamList = [5 10 20 50 100 200];
FOCUS = 1; FOCUS_BEAM = 10;
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
nSat = size(sats,1); nB = numel(beamList); nA = numel(alphaList);
as = pi/180/3600;

%% --------------------------- RUN ---------------------------------------
fprintf('\n=== Points acquisition | N=%d, detElev=%g, PRF=%d, alpha={%.3f,%.3f} ===\n', ...
        Npts, detElevDeg, hw.PRF, alphaList(1), alphaList(2));

L  = chol(diag([sigmaRIC.R^2, sigmaRIC.I^2, sigmaRIC.C^2]), 'lower');
dR = (L * randn(3, Npts))';                 % Npts x 3

acqAll = nan(nSat,nB,nA); rimAll = nan(nSat,nB,nA);
names  = cell(nSat,1);  focusSnap = {};
for s = 1:nSat
    names{s} = sats{s,1};
    cfgBase.mode='tle'; cfgBase.tle = parseTLE(sats{s,2}, sats{s,3});
    pass = findBestPass(cfgBase, station, passOpt);
    cfg = cfgBase; cfg.tSec = pass.tSec; trk = generateTrack(cfg);
    cov = generateCovariance(trk, sigmaRIC);
    fou = projectFoU(cov.SigmaEci, trk.rEci, trk.jd, station, projOpt);

    win = find(fou.altDeg >= detElevDeg); i0 = win(1); i1 = win(end);
    idx = i0:i1; M = numel(idx);
    tGrid  = trk.tSec(idx); meanEl = deg2rad(fou.altDeg(idx)); meanAz = deg2rad(fou.azDeg(idx));
    [~, ipk] = max(meanEl);

    dCross = zeros(Npts,M,'single'); dEl = zeros(Npts,M,'single');
    for jj = 1:M
        m = idx(jj);
        Q = [cov.Rhat(:,m), cov.Ihat(:,m), cov.Chat(:,m)];
        candEci = trk.rEci(:,m) + Q*dR.';
        [caz, cel] = eciToAzElVec(candEci, gmstRad(trk.jd(m)), station);
        dA = mod(caz - meanAz(jj) + pi, 2*pi) - pi;
        dCross(:,jj) = single((dA .* cos(meanEl(jj))).');
        dEl(:,jj)    = single((cel - meanEl(jj)).');
    end
    rho = sqrt(dCross.^2 + dEl.^2);
    RmaxCover = 1.15 * double(max(rho(:)));

    % rim = outer rimFrac by culmination radius
    rc = double(rho(:,ipk)); rcs = sort(rc);
    thr = rcs(max(1, ceil((1-rimFrac)*numel(rcs)))); rimIdx = find(rc >= thr);

    fprintf('\n%-16s window %.1f min (elev %.0f->%.0f->%.0f), a_max %.0f urad, N_rim=%d\n', ...
        sats{s,1}, (tGrid(end)-tGrid(1))/60, rad2deg(meanEl(1)), rad2deg(meanEl(ipk)), ...
        rad2deg(meanEl(end)), max(fou.aRad(idx))*1e6, numel(rimIdx));
    fprintf('   %6s | %-20s | %-20s\n','beam"', ...
        sprintf('alpha=%.3f (Rs=%.2fx)',alphaList(1),2*alphaList(1)), ...
        sprintf('alpha=%.3f (Rs=%.2fx)',alphaList(2),2*alphaList(2)));
    fprintf('   %6s | %9s %9s | %9s %9s\n','','acq','rim','acq','rim');
    for bi = 1:nB
        Rbeam = 0.5*beamList(bi)*as;
        for ai = 1:nA
            sim = spiralAcqSim(tGrid, meanEl, dCross, dEl, rho, Rbeam, alphaList(ai), hw, RmaxCover);
            acqAll(s,bi,ai) = sim.acqRate;
            rimAll(s,bi,ai) = mean(sim.acquired(rimIdx));
            if s==FOCUS && beamList(bi)==FOCUS_BEAM
                focusSnap{ai} = struct('acq',sim.acquired,'sx',sim.sx,'sy',sim.sy, ...
                    'dCrossPk',double(dCross(:,ipk)),'dElPk',double(dEl(:,ipk)), ...
                    'alpha',alphaList(ai),'acqRate',sim.acqRate);
            end
        end
        fprintf('   %6d | %8.1f%% %8.1f%% | %8.1f%% %8.1f%%\n', beamList(bi), ...
            100*acqAll(s,bi,1), 100*rimAll(s,bi,1), 100*acqAll(s,bi,2), 100*rimAll(s,bi,2));
    end
end

%% --------------------------- PLOTS -------------------------------------
mk = {'-o','-s'}; acol = {[0.20 0.45 0.85],[0.85 0.25 0.10]};
alab = arrayfun(@(a) sprintf('alpha=%.3f',a), alphaList, 'UniformOutput',false);

% (1) overall + rim acquisition vs beam, per satellite (2 rows x 3 cols)
fig = figure('Position',[60 60 1280 720],'Color','w','Visible','off');
for s = 1:nSat
    subplot(2,nSat,s); hold on; grid on; box on;
    for ai = 1:nA
        plot(beamList, 100*squeeze(acqAll(s,:,ai)), mk{ai}, 'Color',acol{ai}, ...
             'LineWidth',1.8,'MarkerFaceColor',acol{ai});
    end
    set(gca,'XScale','log'); ylim([0 103]);
    title(sprintf('%s: overall acq', names{s})); xlabel('beam width [arcsec]'); ylabel('acq [%]');
    if s==1, legend(alab,'Location','southeast'); legend boxoff; end

    subplot(2,nSat,nSat+s); hold on; grid on; box on;
    for ai = 1:nA
        plot(beamList, 100*squeeze(rimAll(s,:,ai)), mk{ai}, 'Color',acol{ai}, ...
             'LineWidth',1.8,'MarkerFaceColor',acol{ai});
    end
    set(gca,'XScale','log'); ylim([0 103]);
    title(sprintf('%s: rim acq (outer %.0f%%)', names{s}, 100*rimFrac));
    xlabel('beam width [arcsec]'); ylabel('rim acq [%]');
end
print(fig, fullfile(OUTDIR,'acq_alpha_vs_beam.png'), '-dpng','-r110'); close(fig);

% (2) sky snapshot at culmination for the focus case, both alphas
if numel(focusSnap) == nA
    urad = 1e6;
    fig = figure('Position',[60 60 1200 560],'Color','w','Visible','off');
    for ai = 1:nA
        S = focusSnap{ai};
        subplot(1,nA,ai); hold on; grid on; box on; axis equal;
        plot(S.sx*urad, S.sy*urad, '.', 'Color',[0.75 0.75 0.75], 'MarkerSize',2);
        ac = S.acq; mi = ~S.acq;
        plot(S.dCrossPk(ac)*urad, S.dElPk(ac)*urad, '.', 'Color',[0.20 0.55 0.25],'MarkerSize',3);
        plot(S.dCrossPk(mi)*urad, S.dElPk(mi)*urad, '.', 'Color',[0.85 0.15 0.10],'MarkerSize',4);
        xlabel('cross-elevation [urad]'); ylabel('elevation [urad]');
        title(sprintf('%s  %d arcsec  alpha=%.3f  acq %.1f%%', names{FOCUS}, FOCUS_BEAM, S.alpha, 100*S.acqRate));
    end
    print(fig, fullfile(OUTDIR,'acq_alpha_snapshot.png'), '-dpng','-r110'); close(fig);
end

fprintf('\nFigures: figures/acq_alpha_vs_beam.png, figures/acq_alpha_snapshot.png\n');
fprintf('Done.\n');
