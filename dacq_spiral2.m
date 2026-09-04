function [P, tTrav, info] = dacq_spiral2(pre, cfg, opt)
%==========================================================================
% dacq_spiral2  Archimedean baseline with an EXPLICIT actuator/dwell policy,
%               so the spiral and the serpentine can be held to the same
%               limits and the two kinds of limit can be separated.
%
%   opt.accPolicy : 'legacy'  0.90*a_avail, worst-case cos(el)   (original)
%                   'match'   a_avail, worst-case cos(el)        (= serpentine budget)
%                   'inst'    a_avail, instantaneous cos(el(t))  (fairest)
%   opt.velFac    : fraction of (w_max - centreline rate) usable (legacy 0.80)
%   opt.useDwell  : apply the dwell-spacing cap v <= 2*r_beam*f  (true/false)
%
%   The two caps are physically different things:
%     ACTUATOR  v <= w_max and v <= sqrt(a*R_curv): what the MOUNT can do.
%     DWELL     v <= 2*r_beam*f_dwell: what the SHOT RATE can tile without
%               leaving gaps between consecutive footprints.
%==========================================================================
    if nargin < 3, opt = struct(); end
    if ~isfield(opt,'accPolicy'), opt.accPolicy = 'legacy'; end
    if ~isfield(opt,'velFac'),    opt.velFac    = 0.80;     end
    if ~isfield(opt,'useDwell'),  opt.useDwell  = true;     end
    if ~isfield(opt,'profile'),   opt.profile   = 'naive';  end

    b     = pre.b;
    thMax = pre.Om_max / b;
    th    = linspace(0, thMax, 20000).';
    Rc    = b*(1+th.^2).^1.5 ./ (th.^2 + 2);
    dsdth = b*sqrt(1+th.^2);

    gamAz = max(abs(diff(pre.GAM(:,1))))/pre.dt;
    vSlew = opt.velFac * max(cfg.wMax - gamAz, 0.1*cfg.wMax);
    vDwell= 2*cfg.rBeam*cfg.fDwell;

    switch opt.accPolicy
    case 'legacy', aBase = 0.90*min(pre.a_avail);  cfun = @(t) cos(max(pre.GAM(:,2)));
    case 'match',  aBase =      min(pre.a_avail);  cfun = @(t) cos(max(pre.GAM(:,2)));
    case 'inst',   aBase =      min(pre.a_avail);  cfun = @(t) cos(interp1(pre.tD, pre.GAM(:,2), ...
                                                          min(max(t,0),pre.tD(end)), 'linear'));
    end

    s_arc = cumtrapz(th, dsdth);                     % arclength
    cs = cfun(0)*ones(size(th));                     % first pass: constant
    for it = 1:6                                     % converge cos(el(t))
        v = min(vSlew.*cs, sqrt(aBase*cs.*Rc));
        if opt.useDwell, v = min(v, vDwell); end
        if strcmp(opt.profile,'proper')
            v   = velProfile(v, s_arc, Rc, aBase*cs, pre.dt);
            tth = [0; cumsum(2*diff(s_arc) ./ (v(1:end-1)+v(2:end)))];
        else
            tth = cumtrapz(th, dsdth ./ v);
        end
        if strcmp(opt.accPolicy,'inst')
            csNew = cfun(tth);
            if max(abs(csNew-cs)) < 1e-6, cs = csNew; break; end
            cs = csNew;
        else
            break;
        end
    end

    if strcmp(opt.profile,'proper')
        %  March the traverse forward IN TIME, integrating theta itself.
        %  Two earlier attempts failed and both failures are instructive:
        %  inverting a theta-indexed time table and interpolating linearly is
        %  not self-consistent at the dwell rate (near the end one theta node
        %  spans tens of dwell samples, so the speed is held constant inside
        %  the node and kinks at its edge); and marching arclength against a
        %  fixed s-grid quantises theta (near the centre one dwell advances
        %  less than one grid cell, so theta moves in a staircase). Marching
        %  theta continuously avoids both.
        vcap = v;  dth = th(2)-th(1);
        N = pre.N_max;  thk = zeros(N,1);
        thv = 0;  vv = vcap(1);  tTrav = pre.tD(end);  hit = false;
        for k = 1:N
            q  = min(max(thv/dth + 1, 1), numel(th)-1);
            j  = floor(q);  f = q - j;
            Rk = b*(1+thv^2)^1.5 / (thv^2 + 2);
            ak = aBase * (cs(j) + f*(cs(j+1)-cs(j)));
            if thv >= thMax
                %  Past the commanded radius: brake to a halt at the actuator
                %  limit rather than setting the speed to zero in one sample.
                %  Zeroing it is a velocity step, i.e. an infinite
                %  acceleration -- the very artefact this profile exists to
                %  remove. The small overshoot past Om_max is harmless: it is
                %  outside the cloud.
                if ~hit, tTrav = pre.tD(k); hit = true; end
                vc = 0;
            else
                vc = vcap(j) + f*(vcap(j+1)-vcap(j));
            end
            an = vv^2 / Rk;
            at = sqrt(max(ak^2 - an^2, 0));
            if vc < vv, vv = max(vv - at*pre.dt, 0);
            else,       vv = min(vc, vv + at*pre.dt); end
            thv = thv + vv*pre.dt / (b*sqrt(1+thv^2));
            thk(k) = thv;
        end
    else
        tTrav = tth(end);
        thk = interp1(tth, th, pre.tD, 'linear');
        thk(pre.tD >= tTrav) = thMax;  thk(isnan(thk)) = thMax;
    end

    x  = b*thk .* cos(thk);
    y  = b*thk .* sin(thk);
    el = pre.GAM(:,2) + y;
    az = pre.GAM(:,1) + x ./ cos(el);
    P  = [az, el];

    info.vProfile = v;  info.tth = tth;  info.th = th;
    info.vSlew = vSlew*min(cs); info.vDwell = vDwell;
    info.vCentMax = max(sqrt(aBase*cs.*Rc));
    [~, ib] = min([info.vSlew, vDwell*opt.useDwell + 1e9*(~opt.useDwell), info.vCentMax]);
    bn = {'slew','dwell','centripetal'};  info.bind = bn{ib};
end

% =========================================================================
function v = velProfile(vlim, s, R, acap, dtStep)
% Forward-backward time-optimal velocity profile along a path.
%
%   The naive cap v = sqrt(a*R) bounds only the CENTRIPETAL acceleration and
%   silently allows unbounded TANGENTIAL acceleration -- including an
%   instantaneous stop at the end of the traverse, which is a velocity step,
%   i.e. an infinite acceleration the mount cannot produce.
%
%   Here the TOTAL acceleration is bounded: with a_n = v^2/R the tangential
%   budget left is a_t = sqrt(max(acap^2 - a_n^2, 0)), and the profile is
%   swept forward from v = 0 and backward to v = 0 so the scan both spins up
%   and brakes within the actuator limit.
%   The endpoints are NOT set to exactly zero: as v -> 0 the time integral
%   ds/v stretches without bound and theta(t) becomes ill-conditioned to
%   interpolate, which reintroduces spurious acceleration spikes. They are set
%   instead to a_cap*dt, the largest speed the mount can null in a single
%   dwell interval, so parking at the end is itself within the limit.
    n = numel(s);  v = vlim(:);
    vEnd = min(acap(1)*dtStep, min(vlim));
    v(1) = vEnd;  v(n) = vEnd;
    for i = 1:n-1                                    % forward (spin up)
        an = v(i)^2 / R(i);
        at = sqrt(max(acap(i)^2 - an^2, 0));
        v(i+1) = min(v(i+1), sqrt(v(i)^2 + 2*at*(s(i+1)-s(i))));
    end
    for i = n:-1:2                                   % backward (brake)
        an = v(i)^2 / R(i);
        at = sqrt(max(acap(i)^2 - an^2, 0));
        v(i-1) = min(v(i-1), sqrt(v(i)^2 + 2*at*(s(i)-s(i-1))));
    end
    v = max(v, vEnd);
end
