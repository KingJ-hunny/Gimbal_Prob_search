function drawFreqMarks(tv, yv, tScan, fList)
%==========================================================================
% drawFreqMarks  Overlay vertical scan-completion markers on an az/el-vs-time
%                axis: one dashed line + curve marker + label per control
%                frequency, placed at t = tScan(i) = K/fList(i).
%==========================================================================
    yl = [min(yv) max(yv)];
    if diff(yl) < 1e-6, yl = yl + [-1 1]; end
    for i = 1:numel(fList)
        plot([tScan(i) tScan(i)], yl, '--','Color',[0.55 0.55 0.55],'LineWidth',0.8);
        plot(tScan(i), interp1(tv, yv, tScan(i)), 'o', 'Color',[0.9 0.4 0], ...
             'MarkerFaceColor',[1 0.6 0], 'MarkerSize',5);
        yt = yl(1) + (0.93 - 0.10*mod(i,2))*(yl(2)-yl(1));
        text(tScan(i), yt, sprintf('%dHz',fList(i)), 'FontSize',7, ...
             'HorizontalAlignment','center');
    end
end
