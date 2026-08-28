%==========================================================================
% main_points_acq_test.m
%
%   Monte-Carlo acquisition test with MOVING candidate points, using the
%   constant-linear-velocity Archimedean-spiral scan-point allocation:
%
%       r(theta) = theta * R_s/(2*pi),   R_s = 2*alpha*Rbeam,
%       theta_k  = sqrt(4*pi*k)      =>   r_k = R_s*sqrt(k/pi),
%
%   alpha is the overlap factor (alpha = 1/sqrt(2): full coverage; alpha = 1:
%   footprints tangent, gaps).
%
%   Two candidate distributions over the SAME 3-sigma RIC region are compared:
%     'gaussian' : N(0, Sigma_RIC)                       (probability-weighted)
%     'uniform'  : uniform inside the 3-sigma RIC ellipsoid
%   For the uniform set every hypothesis is equally weighted, so the acquisition
%   rate is the COVERED FRACTION of the 3-sigma region -- a pure coverage/breadth
%   measure, without the Gaussian centre-weighting that flatters the overall
%   number.
%
%   Each candidate is a fixed RIC offset propagated through ECI->ENU and is
%   acquired the first time a scan point lands within Rbeam of its CURRENT
%   position at that point's visit time.
%
%   Swept here: N = 40000 points; distributions {gaussian, uniform};
%   alpha {1/sqrt(2), 1}; all beam widths; detectable elevation {10, 20} deg;
%   3 satellites.
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

Npts        = 40000;
alphaList   = [1/sqrt(2), 1.0];        % overlap factor (R_s = 2*alpha*Rbeam)
distList    = {'gaussian','uniform'};  % candidate distribution over the 3-sigma region
detElevList = [10, 20];                % detectable elevation to start scanning [deg]
kSig        = 3;                        % uniform fills the k-sigma RIC ellipsoid

hw.Vaz = 20; hw.Vel = 10; hw.Aaz = 5; hw.Ael = 2; hw.PRF = 2000;
beamList = [5 10 20 50 100 200];
FOCUS = 1; FOCUS_BEAM = 10; FOCUS_ALPHA = 1; FOCUS_DET = 20;   % for the snapshot
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
nD = numel(distList); nE = numel(detElevList);
as = pi/180/3600;

%% --------------------------- SAMPLE CLOUDS -----------------------------
sig3 = kSig*[sigmaRIC.R, sigmaRIC.I, sigmaRIC.C];
L    = chol(diag([sigmaRIC.R^2, sigmaRIC.I^2, sigmaRIC.C^2]), 'lower');
u    = randn(3, Npts); u = u ./ sqrt(sum(u.^2,1));   % random directions
rr   = rand(1, Npts).^(1/3);                         % radius for uniform-in-ball
dRset = { (L*randn(3,Npts)).', (diag(sig3) * (u.*rr)).' };  % Npts x 3, RIC [km]

%% --------------------------- RUN ---------------------------------------
fprintf('\n=== Points acquisition | N=%d, PRF=%d, alpha={%.3f,%.3f}, detElev={%g,%g} ===\n', ...
        Npts, hw.PRF, alphaList(1), alphaList(2), detElevList(1), detElevList(2));

acqAll = nan(nSat,nB,nA,nD,nE);
names  = cell(nSat,1);  focusSnap = cell(1,nD);
baseDet = min(detElevList);
for s = 1:nSat
    names{s} = sats{s,1};
    cfgBase.mode='tle'; cfgBase.tle = parseTLE(sats{s,2}, sats{s,3});
    pass = findBestPass(cfgBase, station, passOpt);
    cfg = cfgBase; cfg.tSec = pass.tSec; trk = generateTrack(cfg);
    cov = generateCovariance(trk, sigmaRIC);
    fou = projectFoU(cov.SigmaEci, trk.rEci, trk.jd, station, projOpt);

    win = find(fou.altDeg >= baseDet); i0 = win(1); i1 = win(end);
    idx = i0:i1; M = numel(idx);
    tGridB = trk.tSec(idx); meanElB = deg2rad(fou.altDeg(idx)); meanAzB = deg2rad(fou.azDeg(idx));
    [~, ipk] = max(meanElB);

    fprintf('\n%-16s pass window %.1f min (elev %.0f->%.0f->%.0f), a_max %.0f urad\n', ...
        sats{s,1}, (tGridB(end)-tGridB(1))/60, rad2deg(meanElB(1)), rad2deg(meanElB(ipk)), ...
        rad2deg(meanElB(end)), max(fou.aRad(idx))*1e6);

    for di = 1:nD
        dR = dRset{di};
        dCrossB = zeros(Npts,M,'single'); dElB = zeros(Npts,M,'single');
        for jj = 1:M
            m = idx(jj);
            Q = [cov.Rhat(:,m), cov.Ihat(:,m), cov.Chat(:,m)];
            candEci = trk.rEci(:,m) + Q*dR.';
            [caz, cel] = eciToAzElVec(candEci, gmstRad(trk.jd(m)), station);
            dA = mod(caz - meanAzB(jj) + pi, 2*pi) - pi;
            dCrossB(:,jj) = single((dA .* cos(meanElB(jj))).');
            dElB(:,jj)    = single((cel - meanElB(jj)).');
        end
        rhoB = sqrt(dCrossB.^2 + dElB.^2);

        for ei = 1:nE
            sub = find(rad2deg(meanElB) >= detElevList(ei));    % contiguous sub-window
            tGrid = tGridB(sub); meanEl = meanElB(sub);
            dCross = dCrossB(:,sub); dEl = dElB(:,sub); rho = rhoB(:,sub);
            RmaxCover = 1.05 * double(max(rho(:)));
            for bi = 1:nB
                Rbeam = 0.5*beamList(bi)*as;
                for ai = 1:nA
                    sim = spiralAcqSim(tGrid, meanEl, dCross, dEl, rho, Rbeam, alphaList(ai), hw, RmaxCover);
                    acqAll(s,bi,ai,di,ei) = sim.acqRate;
                    if s==FOCUS && beamList(bi)==FOCUS_BEAM && ai==FOCUS_ALPHA && detElevList(ei)==FOCUS_DET
                        focusSnap{di} = struct('acq',sim.acquired,'sx',sim.sx,'sy',sim.sy, ...
                            'dCrossPk',double(dCrossB(:,ipk)),'dElPk',double(dElB(:,ipk)), ...
                            'dist',distList{di},'acqRate',sim.acqRate);
                    end
                end
            end
        end
        % table (this distribution): acq for the two detElev at alpha=1/sqrt(2)
        fprintf('   [%-8s] acq at alpha=%.3f:   %6s | %10s %10s\n', distList{di}, alphaList(1), ...
            'beam"', sprintf('det%g',detElevList(1)), sprintf('det%g',detElevList(2)));
        for bi = 1:nB
            fprintf('              %26s %6d | %9.1f%% %9.1f%%\n', '', beamList(bi), ...
                100*acqAll(s,bi,1,di,1), 100*acqAll(s,bi,1,di,2));
        end
    end
end

%% --------------------------- PLOTS -------------------------------------
% (1) acq vs beam: rows=distribution, cols=satellite; color=alpha, linestyle=detElev
acol = {[0.20 0.45 0.85],[0.85 0.25 0.10]};
els  = {'-','--'};  emk = {'o','s'};
fig = figure('Position',[50 50 1320 760],'Color','w','Visible','off');
for di = 1:nD
  for s = 1:nSat
    subplot(nD,nSat,(di-1)*nSat + s); hold on; grid on; box on;
    for ai = 1:nA
      for ei = 1:nE
        plot(beamList, 100*squeeze(acqAll(s,:,ai,di,ei)), [els{ei} emk{ai}], ...
             'Color',acol{ai}, 'LineWidth',1.6, 'MarkerFaceColor',acol{ai}, 'MarkerSize',4);
      end
    end
    set(gca,'XScale','log'); ylim([0 103]);
    if di==2, ttl = sprintf('%s: uniform (=coverage %%)', names{s});
    else,     ttl = sprintf('%s: gaussian', names{s}); end
    title(ttl); xlabel('beam width [arcsec]'); ylabel('acq [%]');
    if di==1 && s==1
        legend({'a=0.707 det10','a=0.707 det20','a=1 det10','a=1 det20'}, ...
               'Location','southwest','FontSize',7); legend boxoff;
    end
  end
end
print(fig, fullfile(OUTDIR,'acq_dist_det_vs_beam.png'), '-dpng','-r110'); close(fig);

% (2) culmination snapshot (focus case): gaussian vs uniform
if all(~cellfun(@isempty, focusSnap))
    urad = 1e6;
    fig = figure('Position',[60 60 1200 560],'Color','w','Visible','off');
    for di = 1:nD
        S = focusSnap{di};
        subplot(1,nD,di); hold on; grid on; box on; axis equal;
        plot(S.sx*urad, S.sy*urad, '.', 'Color',[0.75 0.75 0.75], 'MarkerSize',2);
        ac = S.acq; mi = ~S.acq;
        plot(S.dCrossPk(ac)*urad, S.dElPk(ac)*urad, '.', 'Color',[0.20 0.55 0.25],'MarkerSize',3);
        plot(S.dCrossPk(mi)*urad, S.dElPk(mi)*urad, '.', 'Color',[0.85 0.15 0.10],'MarkerSize',4);
        xlabel('cross-elevation [urad]'); ylabel('elevation [urad]');
        title(sprintf('%s  %d arcsec  a=%.3f  det%g  %s: acq %.1f%%', ...
            names{FOCUS}, FOCUS_BEAM, alphaList(FOCUS_ALPHA), FOCUS_DET, S.dist, 100*S.acqRate));
    end
    print(fig, fullfile(OUTDIR,'acq_dist_snapshot.png'), '-dpng','-r110'); close(fig);
end

fprintf('\nFigures: figures/acq_dist_det_vs_beam.png, figures/acq_dist_snapshot.png\n');
fprintf('Done.\n');
