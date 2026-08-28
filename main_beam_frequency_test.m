%==========================================================================
% main_beam_frequency_test.m
%
%   Why does the acquisition error grow as the beam width goes 5 -> 200 arcsec,
%   and how does the pulse/dwell frequency (PRF) change that? This isolates the
%   PRF by fixing everything else: alpha = 1/sqrt(2) (full-coverage spacing), a
%   single detectable elevation, and the uniform 3-sigma cloud (so acquisition
%   rate = covered fraction of the 3-sigma region).
%
%   PRF enters only through the scan speed cap: the beam may advance at most one
%   point spacing R_s per pulse, so the on-sky speed is capped at v <= R_s*PRF.
%   A higher PRF therefore lets the spiral sweep FASTER (more of the region in
%   the limited pass time) but does NOT change the spatial spacing R_s. This
%   decomposes the beam-width error into two regimes:
%     * fine beam  -> tiny R_s -> slow, TIME-limited scan  -> PRF helps a lot
%     * coarse beam-> large R_s -> few, sparse dwell points -> PRF cannot help
%
%   Sweeps PRF in {100, 500, 1000, 2000, 3000} Hz over all beam widths, 3 sats.
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
sigmaRIC.R = 0.205; sigmaRIC.I = 0.808; sigmaRIC.C = 0.168;

passOpt.searchHours = 12; passOpt.coarseDtSec = 15; passOpt.fineDtSec = 2; passOpt.maskDeg = 10;
projOpt.nSigma = 3; projOpt.sigmaIsoRad = 0; projOpt.nEllipse = 12;

detElevDeg = 20;              % single detectable elevation (control frequency only)
Npts       = 40000;
alpha      = 1/sqrt(2);       % fixed overlap factor (R_s = sqrt(2)*Rbeam)
kSig       = 3;               % uniform in the 3-sigma RIC ellipsoid

hw.Vaz = 20; hw.Vel = 10; hw.Aaz = 5; hw.Ael = 2;   % PRF set in the loop
prfList  = [100 500 1000 2000 3000];                 % Hz
beamList = [5 10 20 50 100 200];                     % arcsec
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
nSat = size(sats,1); nB = numel(beamList); nP = numel(prfList);
as = pi/180/3600;

% uniform 3-sigma RIC cloud, shared across satellites
sig3 = kSig*[sigmaRIC.R, sigmaRIC.I, sigmaRIC.C];
u  = randn(3, Npts); u = u ./ sqrt(sum(u.^2,1));
rr = rand(1, Npts).^(1/3);
dR = (diag(sig3) * (u.*rr)).';           % Npts x 3

%% --------------------------- RUN ---------------------------------------
fprintf('\n=== Beam vs PRF | N=%d(uniform), detElev=%g, alpha=%.3f, PRF sweep ===\n', ...
        Npts, detElevDeg, alpha);

acqAll = nan(nSat,nB,nP);  names = cell(nSat,1);
for s = 1:nSat
    names{s} = sats{s,1};
    cfgBase.mode='tle'; cfgBase.tle = parseTLE(sats{s,2}, sats{s,3});
    pass = findBestPass(cfgBase, station, passOpt);
    cfg = cfgBase; cfg.tSec = pass.tSec; trk = generateTrack(cfg);
    cov = generateCovariance(trk, sigmaRIC);
    fou = projectFoU(cov.SigmaEci, trk.rEci, trk.jd, station, projOpt);

    win = find(fou.altDeg >= detElevDeg); idx = win(1):win(end); M = numel(idx);
    tGrid = trk.tSec(idx); meanEl = deg2rad(fou.altDeg(idx)); meanAz = deg2rad(fou.azDeg(idx));

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
    rho = sqrt(dCross.^2 + dEl.^2); RmaxCover = 1.05 * double(max(rho(:)));

    fprintf('\n%-16s window %.1f min, a_max %.0f urad\n', ...
        sats{s,1}, (tGrid(end)-tGrid(1))/60, max(fou.aRad(idx))*1e6);
    fprintf('   %6s |', 'beam"'); fprintf(' %7dHz', prfList); fprintf('   (coverage %%)\n');
    for bi = 1:nB
        Rbeam = 0.5*beamList(bi)*as;
        fprintf('   %6d |', beamList(bi));
        for pi_ = 1:nP
            hw.PRF = prfList(pi_);
            sim = spiralAcqSim(tGrid, meanEl, dCross, dEl, rho, Rbeam, alpha, hw, RmaxCover);
            acqAll(s,bi,pi_) = sim.acqRate;
            fprintf(' %8.1f', 100*sim.acqRate);
        end
        fprintf('\n');
    end
end

%% --------------------------- PLOTS -------------------------------------
% (1) coverage vs beam width, one curve per PRF, per satellite
cmap = jet(nP); plab = arrayfun(@(f) sprintf('%d Hz',f), prfList, 'UniformOutput',false);
fig = figure('Position',[60 60 1320 430],'Color','w','Visible','off');
for s = 1:nSat
    subplot(1,nSat,s); hold on; grid on; box on;
    for pi_ = 1:nP
        plot(beamList, 100*squeeze(acqAll(s,:,pi_)), '-o', 'Color',cmap(pi_,:), ...
             'LineWidth',1.6,'MarkerFaceColor',cmap(pi_,:),'MarkerSize',4);
    end
    set(gca,'XScale','log'); ylim([0 103]);
    title(sprintf('%s (uniform, alpha=%.3f)', names{s}, alpha));
    xlabel('beam width [arcsec]'); ylabel('coverage [%]');
    if s==1, legend(plab,'Location','southwest'); legend boxoff; end
end
print(fig, fullfile(OUTDIR,'beam_freq_coverage.png'), '-dpng','-r110'); close(fig);

% (2) coverage vs PRF, one curve per beam width (focus satellite)
sref = 1; cmap2 = jet(nB); blab = arrayfun(@(b) sprintf('%d arcsec',b), beamList, 'UniformOutput',false);
fig = figure('Position',[60 60 720 540],'Color','w','Visible','off');
hold on; grid on; box on;
for bi = 1:nB
    plot(prfList, 100*squeeze(acqAll(sref,bi,:)), '-o', 'Color',cmap2(bi,:), ...
         'LineWidth',1.7,'MarkerFaceColor',cmap2(bi,:));
end
set(gca,'XScale','log'); ylim([0 103]);
xlabel('PRF [Hz]'); ylabel('coverage [%]');
title(sprintf('%s: coverage vs PRF (uniform, alpha=%.3f)', names{sref}, alpha));
legend(blab,'Location','southeast'); legend boxoff;
print(fig, fullfile(OUTDIR,'beam_freq_vs_prf.png'), '-dpng','-r110'); close(fig);

fprintf('\nFigures: figures/beam_freq_coverage.png, figures/beam_freq_vs_prf.png\n');
fprintf('Done.\n');
