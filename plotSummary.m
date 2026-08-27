function plotSummary(R, outdir)
%==========================================================================
% plotSummary  Cross-target comparison: axis ratio vs elevation and the
%              FOU major axis vs range (the 1/range law) for all passes.
%==========================================================================
    urad = 1e6;
    cols = lines(numel(R));

    fig = figure('Position',[100 100 1180 460],'Color','w','Visible','off');

    % ---- axis ratio vs elevation --------------------------------------
    subplot(1,3,1); hold on; grid on; box on;
    for s = 1:numel(R)
        if isempty(R{s}) || ~R{s}.found, continue; end
        plot(R{s}.fou.altDeg, R{s}.fou.k, '-', 'Color',cols(s,:), 'LineWidth',1.6);
    end
    xlabel('elevation [deg]'); ylabel('axis ratio k = b/a');
    title('FOU anisotropy vs elevation');

    % ---- major axis vs range (1/range law) ----------------------------
    subplot(1,3,2); hold on; grid on; box on;
    for s = 1:numel(R)
        if isempty(R{s}) || ~R{s}.found, continue; end
        plot(R{s}.fou.rangeKm, R{s}.fou.aRad*urad, '-', 'Color',cols(s,:), 'LineWidth',1.6);
    end
    xlabel('slant range [km]'); ylabel('major axis a (3-sigma) [urad]');
    title('1/range law');

    % ---- major axis vs elevation --------------------------------------
    subplot(1,3,3); hold on; grid on; box on;
    labs = {};
    for s = 1:numel(R)
        if isempty(R{s}) || ~R{s}.found, continue; end
        plot(R{s}.fou.altDeg, R{s}.fou.aRad*urad, '-', 'Color',cols(s,:), 'LineWidth',1.6);
        labs{end+1} = R{s}.name; %#ok<AGROW>
    end
    xlabel('elevation [deg]'); ylabel('major axis a (3-sigma) [urad]');
    title('FOU size vs elevation');
    legend(labs,'Location','best'); legend boxoff;

    print(fig, fullfile(outdir,'fou_summary.png'), '-dpng', '-r120');
    close(fig);
end
