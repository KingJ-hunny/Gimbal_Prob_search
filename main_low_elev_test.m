%==========================================================================
% main_low_elev_test.m
%
%   Does the near-zenith culmination cause the acquisition failure? Compare a
%   LOW-elevation pass against a near-zenith pass under identical conditions.
%
%   Two targets, same station (Korea), same covariance, same scan:
%     (1) synthetic LEO tuned to culminate ~40 deg over Korea (Kepler mode:
%         Starlink-like a,i so the FOU physics is identical, only the geometry
%         differs), and
%     (2) Starlink 35978, whose best pass culminates ~84 deg (near zenith).
%
%   A low pass has (a) a WEAKER azimuth keyhole (1/cos40=1.3 vs 1/cos84=9.6)
%   and (b) a SMALLER angular FOU (larger range, 1/range), so the fine-beam,
%   time/keyhole-limited misses should shrink; the coarse-beam / alpha=1
%   spatial-gap losses should be roughly unchanged (elevation-independent).
%
%   alpha = 1/sqrt(2), PRF = 2000, N = 40000, both distributions, all beams.
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

detElevList = [10, 20];      % start-scan elevation: 10 = first visibility, 20 = higher start
Npts   = 40000;
alpha  = 1/sqrt(2);
kSig   = 3;
hw.Vaz = 20; hw.Vel = 10; hw.Aaz = 5; hw.Ael = 2; hw.PRF = 2000;
beamList = [5 10 20 50 100 200];
distList = {'gaussian','uniform'};
OUTDIR = 'figures'; if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

epochJD = 2461279.371438;    % same epoch as the Starlink TLE
% target (1): synthetic orbit designed to culminate ~40 deg over Korea
kep40.mode = 'kepler'; kep40.useJ2 = true;
kep40.elements = struct('aKm',6861.2,'ecc',0.001,'incDeg',43.0, ...
                        'raanDeg',284.5,'argpDeg',0,'M0Deg',0,'epochJD',epochJD);
% target (2): Starlink 35978 (near-zenith best pass ~84 deg)
tle84.mode = 'tle';
tle84.tle  = parseTLE('1 66601U 25270P   26238.87143828  .00021431  00000-0  78601-3 0  9993', ...
                      '2 66601  43.0018 175.0289 0001887 251.0337 109.0311 15.27581436 43570');

targets = { 'synthetic 40deg', kep40; 'Starlink 84deg', tle84 };
nT = size(targets,1); nB = numel(beamList); nD = numel(distList); nE = numel(detElevList);
as = pi/180/3600; baseDet = min(detElevList);

%% --------------------------- SAMPLE CLOUDS -----------------------------
sig3 = kSig*[sigmaRIC.R, sigmaRIC.I, sigmaRIC.C];
L    = chol(diag([sigmaRIC.R^2, sigmaRIC.I^2, sigmaRIC.C^2]), 'lower');
uu   = randn(3, Npts); uu = uu ./ sqrt(sum(uu.^2,1)); rr = rand(1, Npts).^(1/3);
dRset = { (L*randn(3,Npts)).', (diag(sig3) * (uu.*rr)).' };

%% --------------------------- RUN ---------------------------------------
fprintf('\n=== Low vs high elevation | N=%d, alpha=%.3f, PRF=%d, detElev={%g,%g} ===\n', ...
        Npts, alpha, hw.PRF, detElevList(1), detElevList(2));
acqAll = nan(nT,nB,nD,nE); names = cell(nT,1); pkEl = zeros(nT,1); aMax = zeros(nT,1);
for tt = 1:nT
    names{tt} = targets{tt,1}; cfgBase = targets{tt,2};
    pass = findBestPass(cfgBase, station, passOpt);
    cfg = cfgBase; cfg.tSec = pass.tSec; trk = generateTrack(cfg);
    cov = generateCovariance(trk, sigmaRIC);
    fou = projectFoU(cov.SigmaEci, trk.rEci, trk.jd, station, projOpt);

    win = find(fou.altDeg >= baseDet); idx = win(1):win(end); M = numel(idx);
    tGridB = trk.tSec(idx); meanElB = deg2rad(fou.altDeg(idx)); meanAzB = deg2rad(fou.azDeg(idx));
    [~, ipk] = max(meanElB); pkEl(tt) = rad2deg(meanElB(ipk)); aMax(tt) = max(fou.aRad(idx))*1e6;

    fprintf('\n%-22s peak elev %.1f deg, range@peak %.0f km, a_max %.0f urad, keyhole 1/cos=%.1f\n', ...
        names{tt}, pkEl(tt), fou.rangeKm(idx(ipk)), aMax(tt), 1/cos(meanElB(ipk)));
    for di = 1:nD
        dR = dRset{di};
        dCrossB = zeros(Npts,M,'single'); dElB = zeros(Npts,M,'single');
        for jj = 1:M
            m = idx(jj); Q = [cov.Rhat(:,m), cov.Ihat(:,m), cov.Chat(:,m)];
            candEci = trk.rEci(:,m) + Q*dR.';
            [caz, cel] = eciToAzElVec(candEci, gmstRad(trk.jd(m)), station);
            dA = mod(caz - meanAzB(jj) + pi, 2*pi) - pi;
            dCrossB(:,jj) = single((dA .* cos(meanElB(jj))).');
            dElB(:,jj)    = single((cel - meanElB(jj)).');
        end
        rhoB = sqrt(dCrossB.^2 + dElB.^2);
        for ei = 1:nE
            sub = find(rad2deg(meanElB) >= detElevList(ei));
            tGrid = tGridB(sub); meanEl = meanElB(sub);
            dCross = dCrossB(:,sub); dEl = dElB(:,sub); rho = rhoB(:,sub);
            RmaxCover = 1.05*double(max(rho(:)));
            for bi = 1:nB
                Rbeam = 0.5*beamList(bi)*as;
                sim = spiralAcqSim(tGrid, meanEl, dCross, dEl, rho, Rbeam, alpha, hw, RmaxCover);
                acqAll(tt,bi,di,ei) = sim.acqRate;
            end
        end
    end
    for ei = 1:nE
        fprintf('   detElev=%2g  %6s | %-10s %-10s\n', detElevList(ei),'beam"','gaussian','uniform');
        for bi = 1:nB
            fprintf('              %6d | %8.1f%% %8.1f%%\n', beamList(bi), ...
                100*acqAll(tt,bi,1,ei), 100*acqAll(tt,bi,2,ei));
        end
    end
end

%% --------------------------- PLOT --------------------------------------
% uniform coverage; columns = detElev; two target curves each
tcol = {[0.20 0.55 0.30],[0.85 0.20 0.15]};   % 40deg green, 84deg red
tlab = arrayfun(@(t) names{t}, 1:nT, 'UniformOutput',false);
fig = figure('Position',[60 60 1150 470],'Color','w','Visible','off');
for ei = 1:nE
    subplot(1,nE,ei); hold on; grid on; box on;
    for tt = 1:nT
        plot(beamList, 100*squeeze(acqAll(tt,:,2,ei)), '-o', 'Color',tcol{tt}, ...
             'LineWidth',1.9,'MarkerFaceColor',tcol{tt});
    end
    set(gca,'XScale','log'); ylim([0 103]);
    title(sprintf('uniform coverage, start elev = %g deg', detElevList(ei)));
    xlabel('beam width [arcsec]'); ylabel('coverage [%]');
    if ei==1, legend(tlab,'Location','southwest'); legend boxoff; end
end
print(fig, fullfile(OUTDIR,'low_vs_high_elev.png'), '-dpng','-r110'); close(fig);

fprintf('\nFigure: figures/low_vs_high_elev.png\n');
fprintf('Done.\n');
