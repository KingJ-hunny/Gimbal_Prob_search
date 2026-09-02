function dacq_plots(what, pre, cfg, OUTDIR, extra)
%==========================================================================
% dacq_plots  Diagnostics for the debris-acquisition study (spec section 8).
%   what = 'centerline' | 'spiral' | 'pattern' | 'profiles' | 'thist'
%          | 'cumulative' | 'control' | 'sweep'
%==========================================================================
    as = 206264.806;                      % rad -> arcsec
    if nargin < 5, extra = []; end
    fig = figure('Position',[40 40 1100 620],'Color','w','Visible','off');

    switch what
    % ------------------------------------------------------------------
    case 'centerline'
        % sample cloud + centreline + normals, in the local tangent frame
        kk = round(linspace(1, pre.Ks, 6));
        subplot(1,2,1); hold on; grid on; box on; axis equal;
        for q = kk
            plot((pre.AZs(:,q)-pre.GAMs(q,1))*cos(pre.GAMs(q,2))*as, ...
                 (pre.ELs(:,q)-pre.GAMs(q,2))*as, '.', 'MarkerSize',2);
        end
        th = linspace(0,2*pi,80);
        plot(pre.Om_max*as*cos(th), pre.Om_max*as*sin(th), 'k--','LineWidth',1.6);
        xlabel('cross-elevation [arcsec]'); ylabel('elevation [arcsec]');
        title(sprintf('cloud about its median (6 epochs), Om max = %.0f arcsec', pre.Om_max*as));

        subplot(1,2,2); hold on; grid on; box on;
        plot(rad2deg(pre.GAMs(:,1)), rad2deg(pre.GAMs(:,2)), 'b-','LineWidth',2);
        q = round(linspace(1,pre.Ks,14));
        L = 0.25;
        for i = q
            plot(rad2deg(pre.GAMs(i,1)) + [0 L*pre.NHATs(i,1)/cos(pre.GAMs(i,2))], ...
                 rad2deg(pre.GAMs(i,2)) + [0 L*pre.NHATs(i,2)], 'r-','LineWidth',1);
        end
        plot(rad2deg(median(pre.AZs,1)), rad2deg(median(pre.ELs,1)), 'k:','LineWidth',1);
        xlabel('azimuth [deg]'); ylabel('elevation [deg]');
        title('centreline (blue) with normals (red); raw median dotted');
        pname = 'dacq_centerline.png';

    % ------------------------------------------------------------------
    case 'spiral'
        subplot(1,2,1); hold on; grid on; box on; axis equal;
        n = min(pre.N_max, numel(pre.tD));
        dx = (pre.Psp(1:n,1)-pre.GAM(1:n,1)).*cos(pre.GAM(1:n,2))*as;
        dy = (pre.Psp(1:n,2)-pre.GAM(1:n,2))*as;
        plot(dx, dy, '-','Color',[0.15 0.35 0.85],'LineWidth',0.4);
        th = linspace(0,2*pi,80);
        plot(pre.Om_max*as*cos(th), pre.Om_max*as*sin(th), 'k--','LineWidth',1.6);
        xlabel('cross-elevation [arcsec]'); ylabel('elevation [arcsec]');
        title(sprintf('baseline spiral about the centreline (c pitch %g)', cfg.cPitch));

        subplot(1,2,2); hold on; grid on; box on;
        ts = sort(pre.spiral.t_hit);
        plot(ts, (1:numel(ts))/numel(ts)*100, 'b-','LineWidth',2);
        plot([pre.t_conv pre.t_conv], [0 100], 'k--','LineWidth',1.2);
        plot([0 pre.t_max], 100*[pre.frac_ref pre.frac_ref], 'r:','LineWidth',1.2);
        xlim([0 pre.t_max]); ylim([0 102]);
        xlabel('time [s]'); ylabel('harvested [%]');
        title(sprintf('spiral coverage: t ref = %.2f s at %.0f%%, traverse %.2f s', ...
              pre.t_ref, 100*pre.frac_ref, pre.t_conv));
        legend({'spiral','traverse end','reference coverage'},'Location','southeast');
        legend boxoff;
        pname = 'dacq_spiral.png';

    % ------------------------------------------------------------------
    case 'pattern'
        n = pre.N_max;
        subplot(1,2,1); hold on; grid on; box on; axis equal;
        dx = (extra.P(1:n,1)-pre.GAM(1:n,1)).*cos(pre.GAM(1:n,2))*as;
        dy = (extra.P(1:n,2)-pre.GAM(1:n,2))*as;
        plot(dx, dy, '-','Color',[0.85 0.25 0.1],'LineWidth',0.4);
        th = linspace(0,2*pi,80);
        plot(pre.Om_max*as*cos(th), pre.Om_max*as*sin(th), 'k--','LineWidth',1.6);
        xlabel('cross-elevation [arcsec]'); ylabel('elevation [arcsec]');
        title('optimised serpentine, relative to the moving centreline');

        subplot(1,2,2); hold on; grid on; box on;
        plot(rad2deg(pre.Psp(:,1)), rad2deg(pre.Psp(:,2)), '-','Color',[0.3 0.5 0.9],'LineWidth',0.4);
        plot(rad2deg(extra.P(:,1)), rad2deg(extra.P(:,2)), '-','Color',[0.85 0.25 0.1],'LineWidth',0.4);
        xlabel('azimuth [deg]'); ylabel('elevation [deg]');
        title('absolute pointing: spiral (blue) vs serpentine (red)');
        pname = 'dacq_pattern.png';

    % ------------------------------------------------------------------
    case 'profiles'
        t = pre.tD;
        subplot(2,2,1); plot(t, extra.A*as,'LineWidth',1.8); grid on; box on;
        hold on; plot(t, pre.Wd*as,'k--','LineWidth',1.2);
        xlabel('t [s]'); ylabel('A [arcsec]'); title('amplitude A vs distribution width Wd');
        subplot(2,2,2); plot(t, rad2deg(extra.a),'LineWidth',1.8); grid on; box on;
        hold on; plot(t, rad2deg(pre.a_avail),'k--','LineWidth',1.2);
        xlabel('t [s]'); ylabel('a [deg/s^2]'); title('acceleration budget vs a avail');
        subplot(2,2,3); plot(t, extra.OM/pre.dt/(2*pi),'LineWidth',1.8); grid on; box on;
        xlabel('t [s]'); ylabel('f [Hz]'); title('instantaneous serpentine frequency');
        subplot(2,2,4); plot(t, sqrt(extra.A.*extra.a)*as,'LineWidth',1.8); grid on; box on;
        xlabel('t [s]'); ylabel('A*Omega [arcsec/s]'); title('transverse sweep speed');
        pname = 'dacq_profiles.png';

    % ------------------------------------------------------------------
    case 'thist'
        subplot(1,2,1);
        hist(pre.spiral.t_hit(pre.spiral.t_hit < pre.t_max), 40);
        grid on; box on; xlabel('t hit [s]'); ylabel('count');
        title(sprintf('spiral (frac %.3f)', pre.spiral.frac_scan));
        subplot(1,2,2);
        hist(extra.t_hit(extra.t_hit < pre.t_max), 40);
        grid on; box on; xlabel('t hit [s]'); ylabel('count');
        title(sprintf('serpentine (frac %.3f)', extra.frac_scan));
        pname = 'dacq_thist.png';

    % ------------------------------------------------------------------
    case 'cumulative'
        hold on; grid on; box on;
        a1 = sort(pre.spiral.t_hit); a2 = sort(extra.t_hit);
        plot(a1, (1:numel(a1))/numel(a1)*100, '-','Color',[0.15 0.35 0.85],'LineWidth',2.2);
        plot(a2, (1:numel(a2))/numel(a2)*100, '-','Color',[0.85 0.25 0.1],'LineWidth',2.2);
        plot([0 pre.t_max], 100*[pre.frac_ref pre.frac_ref],'k:','LineWidth',1.4);
        xlim([0 pre.t_max]); ylim([0 102]);
        xlabel('time [s]'); ylabel('probability mass harvested [%]');
        title(sprintf('coverage vs time  (reference %.0f%%: spiral %.2f s, serpentine %s s)', ...
              100*pre.frac_ref, pre.t_ref, fmtnum(extra.t_reach)));
        legend({'conventional spiral','optimised serpentine','reference coverage'}, ...
               'Location','southeast'); legend boxoff;
        pname = 'dacq_cumulative.png';

    % ------------------------------------------------------------------
    case 'control'
        hold on; grid on; box on;
        u1 = sqrt(sum(pre.spiral.u_hist.^2,2))/(cfg.aMax*pre.dt);
        u2 = sqrt(sum(extra.u_hist.^2,2))/(cfg.aMax*pre.dt);
        plot((0:numel(u1)-1)*pre.dt, u1, '-','Color',[0.15 0.35 0.85],'LineWidth',0.6);
        plot((0:numel(u2)-1)*pre.dt, u2, '-','Color',[0.85 0.25 0.1],'LineWidth',0.6);
        plot([0 pre.t_max],[1 1],'k--','LineWidth',1.2);
        xlabel('time [s]'); ylabel('|u| / (a max * dt)');
        title(sprintf('control effort: spiral J u = %.3f, serpentine J u = %.3f', ...
              mean(u1.^2), mean(u2.^2)));
        legend({'spiral','serpentine','saturation'},'Location','east'); legend boxoff;
        pname = 'dacq_control.png';

    % ------------------------------------------------------------------
    case 'sweep'
        SW = extra;
        subplot(1,2,1); hold on; grid on; box on;
        plot(SW.f, SW.sp, '-o','Color',[0.85 0.25 0.1],'LineWidth',2,'MarkerFaceColor',[0.85 0.25 0.1]);
        plot(SW.f, SW.spq(1,:), ':','Color',[0.85 0.25 0.1]);
        plot(SW.f, SW.spq(2,:), ':','Color',[0.85 0.25 0.1]);
        plot([min(SW.f) max(SW.f)], [1 1], 'k--','LineWidth',1.2);
        set(gca,'XScale','log'); xlabel('f dwell [Hz]'); ylabel('speedup vs spiral');
        title('serpentine speedup (median, IQR dotted)');
        subplot(1,2,2); hold on; grid on; box on;
        plot(SW.f, SW.tconv, '-s','Color',[0.15 0.35 0.85],'LineWidth',2,'MarkerFaceColor',[0.15 0.35 0.85]);
        set(gca,'XScale','log'); xlabel('f dwell [Hz]'); ylabel('spiral traverse t conv [s]');
        for i = 1:numel(SW.f)
            text(SW.f(i), SW.tconv(i), sprintf('  %s', SW.bn{i}), 'FontSize',7);
        end
        title('conventional traverse time and its bottleneck');
        pname = 'dacq_sweep.png';
    end

    print(fig, fullfile(OUTDIR, pname), '-dpng','-r110');
    close(fig);
end
