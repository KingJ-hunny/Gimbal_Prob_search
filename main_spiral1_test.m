%==========================================================================
% main_spiral1_test.m
%
%   Time-marching spiral-vs-FOU race.
%
%   Acquisition begins at first visibility (elevation = 10 deg, rising) and an
%   Archimedean spiral is started from the centre of the mean-trajectory FOU.
%   The spiral winds OUTWARD at the gimbal/PRF-limited rate while the 3-sigma
%   FOU evolves along the pass. At every epoch the spiral radius r_spiral(t) is
%   compared with the current FOU semi-axes a(t), b(t), and the time at which
%   the spiral crosses the FOU boundary is reported for each beam width.
%
%   Note on direction: from 10 deg on the RISING side the FOU angular size
%   GROWS toward culmination (1/range), then shrinks. So the outward spiral
%   and the FOU meet on the descending side (r_spiral reaches a(t) there),
%   unless a wide beam lets the spiral cover the FOU before culmination.
%
%   Depends on the Stage-1/2 pipeline + spiralExitTime.m.
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end

%% ----------------------------- CONFIG ----------------------------------
station.latDeg = 36.3711; station.lonDeg = 127.3617; station.altKm = 0.10;
station.rEcef  = geodetic2ecef(station.latDeg, station.lonDeg, station.altKm);

sigmaRIC.R = 0.205; sigmaRIC.I = 0.808; sigmaRIC.C = 0.168;   % RIC 1-sigma [km]
projOpt.nSigma = 3; projOpt.sigmaIsoRad = 0; projOpt.nEllipse = 40;

passOpt.searchHours = 12; passOpt.coarseDtSec = 15; passOpt.fineDtSec = 1; passOpt.maskDeg = 10;
startElevDeg = 10;                         % begin the spiral at this elevation (rising)

hw.Vaz = 20; hw.Vel = 10;                  % gimbal rate  [deg/s]
hw.Aaz = 5;  hw.Ael = 2;                   % gimbal accel [deg/s^2]
hw.PRF = 2000;                             % pulse rate [Hz]
hw.overlap = 0.0;                          % beam overlap alpha

beamList = [5 10 20 50 100 200];           % beam width [arcsec]
FOCUS = 1;                                 % which satellite to plot in detail
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

%% --------------------------- RUN ---------------------------------------
fprintf('\n=== Spiral-vs-FOU race | start elev %g deg | PRF=%d Hz ===\n', startElevDeg, hw.PRF);
Rall = cell(nSat,1);
for s = 1:nSat
    cfgBase.mode='tle'; cfgBase.tle = parseTLE(sats{s,2}, sats{s,3});
    pass = findBestPass(cfgBase, station, passOpt);
    cfg = cfgBase; cfg.tSec = pass.tSec; trk = generateTrack(cfg);
    cov = generateCovariance(trk, sigmaRIC);
    fou = projectFoU(cov.SigmaEci, trk.rEci, trk.jd, station, projOpt);

    % start from first visibility (rising)
    i0 = find(fou.altDeg >= startElevDeg, 1, 'first');
    t  = trk.tSec(i0:end);
    el = deg2rad(fou.altDeg(i0:end));
    a  = fou.aRad(i0:end);
    b  = fou.bRad(i0:end);
    [~, ipk] = max(el);

    S = struct('name',sats{s,1},'t',t,'el',el,'a',a,'b',b,'ipk',ipk, ...
               'pass',pass,'beam',beamList);
    S.race = cell(1,numel(beamList));
    S.tPer = nan(1,numel(beamList));  S.tFull = nan(1,numel(beamList));

    fprintf('\n%-16s  peak elev %.1f deg, FOU a: %.0f (start) -> %.0f (peak) -> %.0f urad (end)\n', ...
        sats{s,1}, rad2deg(el(ipk)), a(1)*1e6, a(ipk)*1e6, a(end)*1e6);
    fprintf('   %6s %8s %11s %12s %11s  %s\n', ...
        'beam"','ds"','t_first[s]','t_persist[s]','elev@persist','persistently covered?');
    for bi = 1:numel(beamList)
        hw.beamArcsec = beamList(bi);
        race = spiralExitTime(t, el, a, b, hw);
        S.race{bi} = race; S.tPer(bi) = race.tPer; S.tFull(bi) = race.tFull;
        cov_str = ternary(race.persistCovered, 'yes', 'NO (FOU outgrows spiral to end)');
        fprintf('   %6d %8.1f %11s %12s %11s  %s\n', beamList(bi), race.dsArcsec, ...
            fmtnum(race.tFull), fmtnum(race.tPer), fmtnum(race.elPer), cov_str);
    end
    Rall{s} = S;
end

%% --------------------------- PLOTS -------------------------------------
% (1) radius-vs-time race for the focus satellite
S = Rall{FOCUS}; urad = 1e6; trel = S.t - S.t(1);
fig = figure('Position',[80 80 900 620],'Color','w','Visible','off');
hold on; grid on; box on;
hA = plot(trel, S.a*urad, '-', 'Color',[0.85 0.15 0.10], 'LineWidth',2.6);
hB = plot(trel, S.b*urad, '--','Color',[0.85 0.15 0.10], 'LineWidth',1.4);
cmap = jet(numel(beamList)); labs = {};
if exist('OCTAVE_VERSION','builtin'), hbw = zeros(numel(beamList),1); else, hbw = gobjects(numel(beamList),1); end
yTop = 1.3*max(S.a*urad);                         % zoom to the FOU scale
for bi = 1:numel(beamList)
    r = S.race{bi}.r;
    hbw(bi) = plot(trel, r*urad, '-', 'Color',cmap(bi,:), 'LineWidth',1.5);
    if ~isnan(S.race{bi}.idxPer)                  % persistent-coverage crossing
        kf = S.race{bi}.idxPer;
        plot(trel(kf), min(r(kf)*urad,yTop), 'ko', 'MarkerFaceColor',cmap(bi,:), 'MarkerSize',7);
    end
    labs{end+1} = sprintf('%d arcsec', beamList(bi)); %#ok<AGROW>
end
ylim([0 yTop]);
plot(trel(S.ipk)*[1 1], [0 yTop], ':', 'Color',[0.4 0.4 0.4]);   % culmination marker
xlabel('time since acquisition (elev=10 deg) [s]'); ylabel('angular radius [urad]');
title(sprintf('%s: outward spiral radius vs FOU a(t) (dots = persistent coverage)', S.name));
legend([hA;hB;hbw], [{'FOU a (major)','FOU b (minor)'}, labs], 'Location','northeast'); legend boxoff;
text(trel(S.ipk), yTop*0.96, ' culmination', 'Color',[0.4 0.4 0.4]);
print(fig, fullfile(OUTDIR,'spiral1_race.png'), '-dpng','-r120'); close(fig);

% (2) coverage-complete time vs beam width, all satellites
fig = figure('Position',[80 80 760 560],'Color','w','Visible','off');
hold on; grid on; box on; cols = lines(nSat);
for s = 1:nSat
    plot(beamList, Rall{s}.tPer, '-o', 'Color',cols(s,:), 'LineWidth',1.8, 'MarkerFaceColor',cols(s,:));
end
set(gca,'XScale','log');            % log beam width, linear time (t may be 0)
xlabel('beam width [arcsec]'); ylabel('persistent-coverage time [s]');
title('Time until the spiral disk permanently encloses the FOU');
legend(cellfun(@(r) r.name, Rall, 'UniformOutput',false), 'Location','northeast'); legend boxoff;
print(fig, fullfile(OUTDIR,'spiral1_time_vs_beam.png'), '-dpng','-r120'); close(fig);

fprintf('\nFigures: figures/spiral1_race.png, figures/spiral1_time_vs_beam.png\n');
fprintf('Done.\n');
