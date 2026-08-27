%==========================================================================
% main_conventional_scan.m
%
%   Conventional baseline: centre a circular Archimedean spiral on the mean
%   trajectory and scan outward until the whole FOU is covered. Reports the
%   time to complete the scan under the real gimbal velocity/acceleration
%   limits, the 2 kHz pulse rate, and a selectable beam width.
%
%   The FOU is taken at the pass culmination (peak elevation), where it is
%   largest -> the worst case for a single circular scan. Scan time is also
%   swept along the pass to expose the near-zenith azimuth keyhole.
%
%   Depends on the Stage-1 pipeline (generateTrack/Covariance, projectFoU,
%   findBestPass, parseTLE, ...) and scanArchimedeanSpiral.m.
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end

%% ----------------------------- CONFIG ----------------------------------
% Ground station (Korea) — Daejeon / KAIST
station.latDeg = 36.3711; station.lonDeg = 127.3617; station.altKm = 0.10;
station.rEcef  = geodetic2ecef(station.latDeg, station.lonDeg, station.altKm);

% Orbit-prediction error, 1-sigma RIC [km] (in-track dominated; Kim et al.)
sigmaRIC.R = 0.205; sigmaRIC.I = 0.808; sigmaRIC.C = 0.168;

% FOU projection
projOpt.nSigma = 3; projOpt.sigmaIsoRad = 0; projOpt.nEllipse = 80;

% Pass search
passOpt.searchHours = 12; passOpt.coarseDtSec = 15; passOpt.fineDtSec = 2; passOpt.maskDeg = 10;

% ---- Gimbal / beam hardware spec --------------------------------------
hw.Vaz = 20;   hw.Vel = 10;      % max axis rate  [deg/s]
hw.Aaz = 5;    hw.Ael = 2;       % max axis accel [deg/s^2]  (spec is ">=", i.e. conservative)
hw.PRF = 2000;                   % pulse repetition frequency [Hz]
hw.overlap = 0.0;                % beam overlap factor alpha in [0,1)

% Beam widths to evaluate [arcsec] (footprint diameter)
beamList = [5 10 20 50 100 200];

OUTDIR = 'figures';
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

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
fprintf('\n=== Conventional Archimedean-spiral scan | PRF=%d Hz, alpha=%.1f ===\n', hw.PRF, hw.overlap);
fprintf('gimbal: Vaz=%g Vel=%g deg/s, Aaz=%g Ael=%g deg/s^2\n', hw.Vaz,hw.Vel,hw.Aaz,hw.Ael);

Rall = cell(nSat,1);
for s = 1:nSat
    cfgBase.mode = 'tle'; cfgBase.tle = parseTLE(sats{s,2}, sats{s,3});
    pass = findBestPass(cfgBase, station, passOpt);
    cfg = cfgBase; cfg.tSec = pass.tSec; trk = generateTrack(cfg);
    cov = generateCovariance(trk, sigmaRIC);
    fou = projectFoU(cov.SigmaEci, trk.rEci, trk.jd, station, projOpt);
    [~, ipk] = max(fou.altDeg);

    % tracking rates at culmination (finite diff of the mean track)
    dt   = passOpt.fineDtSec;
    azu  = unwrap(deg2rad(fou.azDeg))*180/pi;
    dAz  = gradient(azu, dt);  dEl = gradient(fou.altDeg, dt);   % deg/s

    Rall{s} = struct('name',sats{s,1},'trk',trk,'fou',fou,'ipk',ipk, ...
                     'pass',pass,'dAzPk',dAz(ipk),'dElPk',dEl(ipk));

    aPk = fou.aRad(ipk); el0 = deg2rad(fou.altDeg(ipk));
    fprintf('\n%-16s peak elev %.1f deg, range %.0f km, FOU a=%.0f urad (3sig), pass %.1f min\n', ...
        sats{s,1}, fou.altDeg(ipk), fou.rangeKm(ipk), aPk*1e6, pass.durationSec/60);
    fprintf('   tracking rate at peak: dAz=%.2f deg/s (limit %g), dEl=%.2f deg/s (limit %g)%s\n', ...
        dAz(ipk), hw.Vaz, dEl(ipk), hw.Vel, ternary(abs(dAz(ipk))>hw.Vaz,'  <-- AZ KEYHOLE > limit',''));
    fprintf('   %6s %8s %7s %9s %9s %10s  %-10s %s\n', ...
        'beam"','ds"','Nturn','T_vel[s]','T_full[s]','T_dwellLB','bind','feasible?');

    Tfull = zeros(size(beamList));
    for bi = 1:numel(beamList)
        hw.beamArcsec = beamList(bi);
        o = scanArchimedeanSpiral(aPk, el0, hw);
        Tfull(bi) = o.T_full;
        [~, wmax] = max([o.bindFrac.v_az o.bindFrac.v_el o.bindFrac.v_sky_PRF o.bindFrac.a_az o.bindFrac.a_el]);
        feas = ternary(o.T_full <= pass.durationSec, 'yes', 'NO (> pass)');
        fprintf('   %6d %8.1f %7.0f %9.2f %9.2f %10.3f  %-10s %s\n', ...
            beamList(bi), o.dsArcsec, o.Nturn, o.T_vel, o.T_full, o.T_dwellLB, o.bindNames{wmax}, feas);
    end
    Rall{s}.Tfull = Tfull;
end

%% --------------------------- PLOTS -------------------------------------
cols = lines(nSat);

% (1) scan time vs beam width, all satellites, with pass-duration references
fig = figure('Position',[80 80 760 560],'Color','w','Visible','off');
hold on; grid on; box on;
for s = 1:nSat
    loglog(beamList, Rall{s}.Tfull, '-o', 'Color',cols(s,:), 'LineWidth',1.8, 'MarkerFaceColor',cols(s,:));
end
set(gca,'XScale','log','YScale','log');
for s = 1:nSat
    plot(xlim, Rall{s}.pass.durationSec*[1 1], '--', 'Color',cols(s,:), 'LineWidth',1.0);
end
xlabel('beam width [arcsec]'); ylabel('scan time to cover FOU [s]');
title('Conventional Archimedean scan time vs beam width (at culmination)');
legend(cellfun(@(r) r.name, Rall, 'UniformOutput',false), 'Location','northeast'); legend boxoff;
text(beamList(1), Rall{1}.pass.durationSec, '  dashed = pass duration', 'VerticalAlignment','bottom');
print(fig, fullfile(OUTDIR,'scan_time_vs_beam.png'), '-dpng','-r120'); close(fig);

% (2) scan time vs elevation along the pass (Starlink 35978), per beam width
sref = 1; r = Rall{sref}; fou = r.fou;
ne = 24; ie = round(linspace(1, numel(fou.altDeg), ne));
fig = figure('Position',[80 80 780 560],'Color','w','Visible','off');
hold on; grid on; box on;
cmap = jet(numel(beamList));
labs = {};
for bi = 1:numel(beamList)
    hw.beamArcsec = beamList(bi);
    Tel = zeros(1,ne);
    for q = 1:ne
        k = ie(q);
        o = scanArchimedeanSpiral(fou.aRad(k), deg2rad(fou.altDeg(k)), hw);
        Tel(q) = o.T_full;
    end
    plot(fou.altDeg(ie), Tel, '-o', 'Color',cmap(bi,:), 'LineWidth',1.5, 'MarkerSize',3);
    labs{end+1} = sprintf('%d arcsec', beamList(bi)); %#ok<AGROW>
end
set(gca,'YScale','log');
xlabel('elevation [deg]'); ylabel('scan time [s] (log)');
title(sprintf('%s: scan time vs elevation (FOU growth + az keyhole)', r.name));
legend(labs,'Location','northwest'); legend boxoff;
print(fig, fullfile(OUTDIR,'scan_time_vs_elev.png'), '-dpng','-r120'); close(fig);

% (3) spiral geometry over the FOU at culmination (beam = 20")
hw.beamArcsec = 20;
r = Rall{1}; fou = r.fou; ipk = r.ipk;
o = scanArchimedeanSpiral(fou.aRad(ipk), deg2rad(fou.altDeg(ipk)), hw);
fig = figure('Position',[80 80 1100 460],'Color','w','Visible','off');
subplot(1,2,1); hold on; grid on; box on; axis equal;
th = linspace(0,2*pi,160);
[V,D] = eig(fou.SigmaAng(:,:,ipk)); lam = max(diag(D),0);
ex = projOpt.nSigma*(V(1,1)*sqrt(lam(1))*cos(th)+V(1,2)*sqrt(lam(2))*sin(th));
ey = projOpt.nSigma*(V(2,1)*sqrt(lam(1))*cos(th)+V(2,2)*sqrt(lam(2))*sin(th));
plot(o.x*1e6, o.y*1e6, '-', 'Color',[0.2 0.45 0.85], 'LineWidth',0.8);
plot(ex*1e6, ey*1e6, '-', 'Color',[0.85 0.2 0.1], 'LineWidth',2.0);
xlabel('cross-elevation [urad]'); ylabel('elevation [urad]');
title(sprintf('%s: spiral (20 arcsec) over FOU, a=%.0f urad', r.name, fou.aRad(ipk)*1e6));
legend({'spiral path','FOU 3-sigma'},'Location','southoutside'); legend boxoff;

subplot(1,2,2); hold on; grid on; box on;
turns = (o.phi - o.phi(1))/(2*pi);
plot(turns, o.MVC,    '--', 'Color',[0.5 0.5 0.5], 'LineWidth',1.0);
plot(turns, o.phidot, '-',  'Color',[0.1 0.1 0.1], 'LineWidth',1.5);
xlabel('spiral turn'); ylabel('d\phi/dt [rad/s]');
title(sprintf('speed profile (T_{full}=%.1f s)', o.T_full));
legend({'MVC (limit)','achieved'},'Location','northeast'); legend boxoff;
print(fig, fullfile(OUTDIR,'scan_spiral_demo.png'), '-dpng','-r120'); close(fig);

fprintf('\nFigures written to ./%s/ (scan_time_vs_beam, scan_time_vs_elev, scan_spiral_demo)\n', OUTDIR);
fprintf('Done.\n');
