function out = spiralAcqSim(tGrid, elMean, dCross, dEl, rho, Rbeam, alpha, hw, RmaxCover)
%==========================================================================
% spiralAcqSim  Monte-Carlo acquisition of MOVING candidate points by an
%               outward, overlapping Archimedean-spiral scan.
%
%   The gimbal tracks the mean trajectory and superimposes an Archimedean
%   spiral of scan points placed by the constant-linear-velocity rule:
%       r(theta) = theta * R_s/(2*pi),   R_s = 2*alpha*Rbeam,
%       theta_k  = sqrt(4*pi*k),  =>  r_k = R_s*sqrt(k/pi).
%   One turn (theta += 2*pi) advances r by exactly R_s, and the along-arc
%   point spacing is also R_s, so R_s is both the arm pitch and the point
%   spacing. alpha is the overlap factor: alpha = 1/sqrt(2) gives R_s =
%   sqrt(2)*Rbeam (full coverage, the diagonal interstice just reaches Rbeam),
%   alpha = 1 gives R_s = 2*Rbeam (footprints tangent, leaving diagonal gaps).
%   The beam visits point k at time t_k, obtained
%   by traversing the spiral under the gimbal/PRF speed limits (velocity,
%   PRF no-gap, and centripetal acceleration, all with the 1/cos(el) azimuth
%   keyhole), exactly as in the earlier stages.
%
%   Candidate points MOVE: their local offset from the mean pointing,
%   (dCross(i,t), dEl(i,t)) with radius rho(i,t), is supplied on the time
%   grid (built in the caller from N covariance samples propagated through
%   ECI->ENU). Candidate i is ACQUIRED the first time a scan point falls
%   within Rbeam of the candidate's CURRENT position at that scan point's
%   visit time -- so a candidate can be missed if the moving/growing cloud is
%   never coincident with the beam when the beam is there (an angular-timing
%   miss the deterministic overlap cannot rescue).
%
%   Given:
%     tGrid   1xM   window times [s] (mean elevation >= detectable elevation)
%     elMean  1xM   mean-pointing elevation [rad]
%     dCross  NxM   candidate cross-elevation offset from mean [rad]
%     dEl     NxM   candidate elevation offset from mean [rad]
%     rho     NxM   candidate radial offset sqrt(dCross^2+dEl^2) [rad]
%     Rbeam   scalar effective beam radius [rad]
%     alpha   spiral overlap/spacing factor (pitch = alpha*Rbeam)
%     hw      .Vaz .Vel [deg/s], .Aaz .Ael [deg/s^2], .PRF [Hz]
%     RmaxCover  radius to scan out to [rad] (e.g. 1.15 * max candidate radius)
%
%   Returned struct out:
%     .acquired   Nx1 logical
%     .acqRate    fraction acquired
%     .tHit       Nx1 acquisition time [s from window start] (NaN if missed)
%     .Kpts       number of scan points actually visited within the window
%     .reached    true if the spiral reached RmaxCover before the window ended
%     .sx .sy .st subsampled scan-point offsets [rad] and times [s] (for plots)
%     .rBeamGrid  1xM spiral radius resampled on tGrid [rad] (for plots)
%==========================================================================

    deg = pi/180;
    Vaz = hw.Vaz*deg; Vel = hw.Vel*deg; Aaz = hw.Aaz*deg; Ael = hw.Ael*deg;
    PRF = hw.PRF;
    [N, M] = size(rho);
    t0 = tGrid(1); tEnd = tGrid(end);
    dtG = tGrid(2) - tGrid(1);

    % ---- spiral scan-point geometry (vectorised) ---------------------
    Rs = 2*alpha*Rbeam;                       % arm pitch = along-arc spacing
    p  = Rs;                                  % (speed cap below uses p as R_s)
    b  = p/(2*pi);                            % r(theta) = (R_s/2pi)*theta
    Kmax = max(4, ceil(pi*(RmaxCover/p)^2));  % r_K ~ RmaxCover
    k  = (1:Kmax);
    th = sqrt(4*pi*k);
    rk = b*th;                                % = p*sqrt(k/pi)
    xk = rk.*cos(th);
    yk = rk.*sin(th);

    % ---- visit time t_k by traversing the spiral (2-pass el estimate) --
    elFun = @(tt) interp1(tGrid, elMean, min(max(tt,t0),tEnd), 'linear');
    tk = t0 + cumsum(repmat(p,1,Kmax) ./ speedCap(rk, elMean(1)*ones(1,Kmax)));
    for it = 1:2
        elk = elFun(tk);
        dtk = p ./ speedCap(rk, elk);
        tk  = t0 + cumsum(dtk);
    end

    reached = (tk(end) <= tEnd);
    keep = tk <= tEnd;                        % scan points reached within window
    xk = xk(keep); yk = yk(keep); rk = rk(keep); tk = tk(keep);
    Kpts = numel(tk);

    % ---- acquisition: interval-grouped ring test (exact angular) -------
    acquired = false(N,1);
    tHit     = nan(N,1);
    if Kpts > 0
        mIdx = min(max(floor((tk - t0)/dtG) + 1, 1), M-1);  % grid interval per point
        [mSort, ord] = sort(mIdx);
        xs = xk(ord); ys = yk(ord); rs = rk(ord); ts = tk(ord);
        edges = [find(diff(mSort)~=0), numel(mSort)];       % last index of each block
        b0 = 1;
        Rb2 = Rbeam^2;
        for e = 1:numel(edges)
            b1 = edges(e); m = mSort(b0);
            xkm = xs(b0:b1); ykm = ys(b0:b1); rkm = rs(b0:b1); tkm = ts(b0:b1);
            b0 = b1 + 1;
            rmin = min(rkm) - Rbeam; rmax = max(rkm) + Rbeam;
            ring = find(rho(:,m) >= rmin & rho(:,m) <= rmax & ~acquired);
            if isempty(ring), continue; end
            cx = dCross(ring, m); cy = dEl(ring, m);         % n_ring x 1
            % distance^2 between every scan point (rows) and ring candidate (cols)
            D2 = (xkm(:) - cx(:).').^2 + (ykm(:) - cy(:).').^2;  % P x n_ring
            hitMask = D2 <= Rb2;
            hitAny = any(hitMask, 1);                        % 1 x n_ring
            if any(hitAny)
                hitIdx = ring(hitAny);
                acquired(hitIdx) = true;
                % first scan point (in time) that hit each newly-acquired cand
                for jj = find(hitAny)
                    firstP = find(hitMask(:,jj), 1, 'first');
                    tHit(ring(jj)) = tkm(firstP) - t0;
                end
            end
        end
    end

    % ---- outputs / plotting aids --------------------------------------
    out.acquired = acquired;
    out.acqRate  = mean(acquired);
    out.tHit     = tHit;
    out.Kpts     = Kpts;
    out.reached  = reached;

    nsub = min(Kpts, 4000);                    % subsample scan points for plots
    if Kpts > 0
        idx = round(linspace(1, Kpts, nsub));
        out.sx = xk(idx); out.sy = yk(idx); out.st = tk(idx);
        out.rBeamGrid = interp1(tk, rk, tGrid, 'linear', NaN);
    else
        out.sx = []; out.sy = []; out.st = []; out.rBeamGrid = nan(1,M);
    end

    % ---- nested speed cap ---------------------------------------------
    function v = speedCap(rr, el)
        ce   = max(cos(el), 1e-3);
        vPRF = p*PRF;
        vVel = min(Vel, Vaz.*ce);
        vAcc = min(sqrt(Ael*max(rr,b)), sqrt(Aaz*max(rr,b).*ce));
        v    = min(min(vPRF, vVel), vAcc);
    end
end
