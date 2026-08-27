function plotSatelliteFOU(res, projOpt, outdir)
%==========================================================================
% plotSatelliteFOU  Six-panel figure summarising the projected FOU for one
%                   satellite pass. See main_projection_test.m.
%==========================================================================
    fou = res.fou; trk = res.trk;
    urad = 1e6;
    tmin = (trk.tSec - trk.tSec(1)) / 60;
    ipk  = res.ipk;

    % representative epochs: rise / peak / set
    idx3   = [1, ipk, numel(fou.altDeg)];
    names3 = {'rise','peak','set'};
    col3   = [0.20 0.45 0.85; 0.85 0.35 0.15; 0.30 0.65 0.30];

    fig = figure('Position',[80 80 1280 760],'Color','w','Visible','off');

    th = linspace(0,2*pi,120);

    % ---- (1) Sky track (pass geometry) --------------------------------
    % Azimuth is unwrapped so the near-zenith fast az sweep reads as a smooth
    % curve rather than a 0/360 jump. (The exaggeration input is unused now:
    % the real FOU is ~tens of urad and is shown to scale in panels 2 & 6.)
    subplot(2,3,1); hold on; grid on; box on;
    azu = unwrap(deg2rad(fou.azDeg)) * 180/pi;
    plot(azu, fou.altDeg, '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1.4);
    for m = 1:3
        plot(azu(idx3(m)), fou.altDeg(idx3(m)), 'o', ...
             'MarkerFaceColor',col3(m,:), 'MarkerEdgeColor','k', 'MarkerSize',7);
    end
    xlabel('azimuth (unwrapped) [deg]'); ylabel('elevation [deg]');
    title('sky track (pass geometry)'); axis tight;

    % ---- (2) FOU ellipses along the pass, true urad, coloured by elev --
    % Overlaying many true-scale ellipses shows both effects at once: they
    % GROW toward high elevation (1/range) and ROTATE as the track curves.
    subplot(2,3,2); hold on; grid on; box on; axis equal;
    nOv = 12;
    io  = round(linspace(1, numel(fou.altDeg), nOv));
    cmap = jet(nOv);
    for jj = 1:nOv
        [bx, by] = ellipseXY(fou.SigmaAng(:,:,io(jj)), projOpt.nSigma, th);
        plot(bx*urad, by*urad, '-', 'Color',cmap(jj,:), 'LineWidth',1.2);
    end
    xlabel('cross-elevation [urad]'); ylabel('elevation [urad]');
    title('FOU vs elevation (blue low -> red high)');
    colormap(jet);
    try
        caxis([min(fou.altDeg) max(fou.altDeg)]);
        hcb = colorbar; ylabel(hcb, 'elevation [deg]');
    catch
    end

    % ---- (3) FOU semi-axes vs elevation -------------------------------
    subplot(2,3,3); hold on; grid on; box on;
    plot(fou.altDeg, fou.aRad*urad, '-', 'Color',[0.85 0.35 0.15],'LineWidth',1.6);
    plot(fou.altDeg, fou.bRad*urad, '-', 'Color',[0.20 0.45 0.85],'LineWidth',1.6);
    xlabel('elevation [deg]'); ylabel('angular half-axis [urad]');
    title('FOU semi-axes (3-sigma)');
    legend({'a (major, in-track)','b (minor, cross)'},'Location','best'); legend boxoff;

    % ---- (4) Axis ratio vs elevation ----------------------------------
    subplot(2,3,4); hold on; grid on; box on;
    plot(fou.altDeg, fou.k, '-', 'Color',[0.35 0.20 0.60],'LineWidth',1.6);
    xlabel('elevation [deg]'); ylabel('axis ratio k = b/a');
    title('FOU anisotropy');

    % ---- (5) The 1/range law and its geometric deviation --------------
    % Compare the fully-projected major axis a against the scalar law
    % 3*sigma_I/range. They coincide near zenith (in-track projects fully
    % onto the sky) and separate at low elevation, where part of the
    % in-track error falls along the line of sight and is lost to pointing.
    subplot(2,3,5); hold on; grid on; box on;
    aPred = projOpt.nSigma * res.cov.sigmaRIC(2) ./ fou.rangeKm;   % 3*sigma_I/range
    plot(fou.altDeg, fou.aRad*urad, '-',  'Color',[0.85 0.35 0.15], 'LineWidth',1.8);
    plot(fou.altDeg, aPred*urad,    '--', 'Color',[0.20 0.30 0.55], 'LineWidth',1.4);
    xlabel('elevation [deg]'); ylabel('major axis a [urad]');
    title('1/range law: projected vs 3*sigmaI/range');
    legend({'a (projected)','3*sigmaI/range'},'Location','northwest'); legend boxoff;

    % ---- (6) Ellipse tilt vs elevation --------------------------------
    subplot(2,3,6); hold on; grid on; box on;
    plot(fou.altDeg, fou.angleDeg, '-', 'Color',[0.15 0.55 0.55],'LineWidth',1.6);
    xlabel('elevation [deg]'); ylabel('major-axis tilt [deg]');
    title('FOU orientation (from cross-El axis)');

    sgtitleCompat(sprintf('%s   |   %s', res.name, res.label));

    fname = fullfile(outdir, sprintf('fou_%s.png', safeName(res.name)));
    print(fig, fname, '-dpng', '-r120');
    close(fig);
end

% -------------------------------------------------------------------------
function [bx, by] = ellipseXY(Sang, nSig, th)
% Boundary of an nSig covariance ellipse from a 2x2 covariance Sang.
    [V, D] = eig(Sang);
    lam = max(diag(D), 0);
    bx = nSig*(V(1,1)*sqrt(lam(1))*cos(th) + V(1,2)*sqrt(lam(2))*sin(th));
    by = nSig*(V(2,1)*sqrt(lam(1))*cos(th) + V(2,2)*sqrt(lam(2))*sin(th));
end

% -------------------------------------------------------------------------
function s = safeName(s)
    s = regexprep(s, '[^A-Za-z0-9]+', '_');
end

function sgtitleCompat(txt)
% sgtitle exists in modern MATLAB and Octave; fall back to an annotation.
    try
        sgtitle(txt, 'FontWeight','bold');
    catch
        try, annotation('textbox',[0 0.95 1 0.05],'String',txt, ...
                'HorizontalAlignment','center','EdgeColor','none','FontWeight','bold'); end
    end
end
