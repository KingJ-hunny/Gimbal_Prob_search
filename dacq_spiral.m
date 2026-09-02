function [P, tTrav] = dacq_spiral(pre, cfg)
%==========================================================================
% dacq_spiral  Conventional Archimedean-spiral baseline, centred on the
%              moving centreline and time-parameterised under the real
%              limits.
%
%   r = b*theta with b = c_pitch*(2 r_beam)/(2 pi), so adjacent turns are one
%   beam diameter apart. The along-path speed is the minimum of
%     * the slew cap        (w_max, less the centreline feedforward, and
%                            de-rated by cos(el) because an on-sky offset
%                            costs 1/cos(el) on the azimuth axis),
%     * the dwell-spacing cap 2*r_beam*f_dwell  (consecutive dwells must
%                            overlap, otherwise the path itself has gaps),
%     * the centripetal cap sqrt(a_avail * R_curv), which is what slows the
%                            spiral near its centre.
%   Whichever cap binds is exactly the "dwell vs slew" bottleneck.
%==========================================================================
    b     = pre.b;
    thMax = pre.Om_max / b;
    th    = linspace(0, thMax, 20000).';
    Rc    = b*(1+th.^2).^1.5 ./ (th.^2 + 2);        % radius of curvature

    cmin  = cos(max(pre.GAM(:,2)));                  % worst-case az inflation
    gamAz = max(abs(diff(pre.GAM(:,1))))/pre.dt;     % centreline az rate
    vSlew = 0.80 * max(cfg.wMax - gamAz, 0.1*cfg.wMax) * cmin;
    vDwell= 2*cfg.rBeam*cfg.fDwell;
    aCap  = 0.90 * min(pre.a_avail) * cmin;

    v     = min(min(vSlew, vDwell), sqrt(aCap*Rc));
    dtdth = b*sqrt(1+th.^2) ./ v;
    tth   = cumtrapz(th, dtdth);

    tTrav = tth(end);                                % one full traverse
    thk = interp1(tth, th, pre.tD, 'linear');
    thk(pre.tD >= tTrav) = thMax;                    % hold once traversed
    thk(isnan(thk)) = thMax;

    x  = b*thk .* cos(thk);
    y  = b*thk .* sin(thk);
    el = pre.GAM(:,2) + y;
    az = pre.GAM(:,1) + x ./ cos(el);
    P  = [az, el];
end
