function pre = dacq_precompute(cfg, pass, nPts, verbose)
%==========================================================================
% dacq_precompute  LAYER 0 of the debris acquisition study.
%
%   A. sample the TLE-error distribution in the RIC frame at epoch
%   B. coarse pass scan -> visibility window, Om_max, nu
%   C. analytic Archimedean-spiral sizing -> N_spiral, bottleneck
%   D. fine propagation of the cloud onto a storage grid (az/el)
%   E. centreline / normal / distribution width / available acceleration
%   F. dwell-grid resampling + Lagrange basis
%   G. spiral baseline simulation -> t_conv, t_max, N_max
%
%   The sample-point az/el are kept on a coarse storage grid (cfg.storeDt)
%   and interpolated to the dwell grid inside dacq_evaluate: the cloud's
%   angular acceleration is ~20 arcsec/s^2, so 0.2 s linear interpolation
%   costs <0.1 arcsec against a 20 arcsec beam, and the memory becomes
%   independent of f_dwell.
%==========================================================================
    if nargin < 3 || isempty(nPts),    nPts = cfg.nPoints; end
    if nargin < 4, verbose = true; end
    dt = 1/cfg.fDwell;

    % =====================================================================
    %  A. sample points from the RIC error distribution
    % =====================================================================
    r0 = pass.r0(:); v0 = pass.v0(:);
    Rh = r0/norm(r0);
    Ch = cross(r0,v0); Ch = Ch/norm(Ch);
    Th = cross(Ch,Rh);
    Mric = [Rh, Th, Ch];

    rand('seed', cfg.seed); randn('seed', cfg.seed);
    dR = (cfg.sigPosRIC(:) * ones(1,nPts)) .* randn(3,nPts);
    dV = (cfg.sigVelRIC(:) * ones(1,nPts)) .* randn(3,nPts);
    sr = r0*ones(1,nPts) + Mric*dR;
    sv = v0*ones(1,nPts) + Mric*dV;

    % =====================================================================
    %  B. coarse scan of the pass -> window and Om_max
    % =====================================================================
    tC0 = max(0, pass.tPeak - 600);
    [rc, vc] = dacq_rk4J2(sr, sv, cfg.rk4Step, round(tC0/cfg.rk4Step));
    tCo = 0:cfg.coarseDt:1200;
    [AZc, ELc] = dacq_trackAzEl(rc, vc, cfg.epochJD + tC0/86400, tCo, ...
                                cfg.station, cfg.rk4Step);

    AZc = unwrapRows(AZc);             % azimuth continuity before any metric

    % the WHOLE cloud must clear the mask, not just its median: a trailing
    % member still below 20 deg cannot be pointed at even in principle
    elMin = min(ELc, [], 1);
    elMed = median(ELc, 1);
    vis   = find(elMin >= deg2rad(cfg.elMaskDeg));
    if isempty(vis), error('dacq_precompute:novis','no visible arc above the mask'); end
    tStart = tC0 + tCo(vis(1));
    tEnd   = tC0 + tCo(vis(end));

    % per-epoch angular radius of the cloud about its median (cos(el) metric)
    radK = zeros(1, numel(tCo));
    for k = vis
        caz = median(AZc(:,k)); cel = median(ELc(:,k));
        radK(k) = max(hypot((AZc(:,k)-caz)*cos(cel), ELc(:,k)-cel));
    end

    % =====================================================================
    %  C. spiral sizing -- Om_max over the SCAN HORIZON, not the whole pass
    % =====================================================================
    %  The cloud is essentially "the target at t +/- dT" with dT set by the
    %  in-track spread, so its sky extent is dT * (sky angular rate). That
    %  rate triples between rise and culmination, so sizing the spiral on the
    %  full visible arc would oversize it ~3x for a scan that only lasts tens
    %  of seconds. Horizon and Om_max depend on each other -> fixed point.
    %  Three caps, not two. Besides the dwell-spacing cap (2*r_beam*f) and
    %  the slew cap (w_max), an Archimedean spiral of radius Om_max cannot be
    %  flown faster than the CENTRIPETAL limit v <= sqrt(a*R_curv): rounding a
    %  0.76 deg circle at 5 deg/s^2 caps v near 1.9 deg/s, well under the
    %  5.6 deg/s the dwell spacing would allow at 500 Hz. Integrating
    %  ds/sqrt(a*b*theta) gives t_cent = sqrt(b/a)*(2/3)*theta_max^1.5.
    b     = cfg.cPitch * (2*cfg.rBeam) / (2*pi);
    T_est = 60;
    for it = 1:10
        kw       = vis(tCo(vis) <= tCo(vis(1)) + T_est);
        Om_max   = 1.02 * max(radK(kw));
        cosEst   = cos(max(elMed(kw)));
        N_spiral = ceil( pi*Om_max^2 / (2*cfg.rBeam*cfg.cPitch)^2 );
        thMax    = Om_max / b;
        L_spiral = 0.5 * b * thMax^2;
        aCapEst  = 0.90 * cfg.rho * cfg.aMax * cosEst;
        t_dwell  = N_spiral / cfg.fDwell;
        t_slew   = L_spiral / cfg.wMax + cfg.wMax / cfg.aMax;
        t_cent   = sqrt(b/aCapEst) * (2/3) * thMax^1.5;
        T_new    = min(cfg.horizonMul * max([t_dwell, t_slew, t_cent]), ...
                       0.90*(tEnd - tStart));
        done     = abs(T_new - T_est) < 0.03*T_est;
        T_est    = T_new;
        if done, break; end
    end
    nu = Om_max / cfg.rBeam;
    [~, ib_] = max([t_dwell, t_slew, t_cent]);
    bn = {'dwell','slew','centripetal'};  bottleneck = bn{ib_};

    if verbose
        fprintf('  visible arc  : %.1f s (el mask %g deg), peak el %.1f deg\n', ...
                tEnd-tStart, cfg.elMaskDeg, pass.elPeakDeg);
        fprintf('  Om_max = %.0f arcsec,  nu = Om_max/r_beam = %.1f\n', ...
                Om_max*206264.806, nu);
        fprintf('  N_spiral = %d   t_dwell = %.1f   t_slew = %.1f   t_cent = %.1f s -> bottleneck = %s\n', ...
                N_spiral, t_dwell, t_slew, t_cent, bottleneck);
    end
    if ~(nu > 5 && nu < 500)
        warning('dacq:nu','nu = %.1f outside (5,500) -- rescale the scenario', nu);
    end

    % =====================================================================
    %  D. fine propagation onto the storage grid
    % =====================================================================
    %  dExt: how many seconds of centreline the cloud spans along-track.
    %  The scan path is extended by this much at both ends (see dacq_setgrid).
    k1 = vis(1);  k2 = min(k1+1, numel(tCo));
    c1 = [median(AZc(:,k1)) median(ELc(:,k1))];
    c2 = [median(AZc(:,k2)) median(ELc(:,k2))];
    skySpd = hypot((c2(1)-c1(1))*cos(c1(2)), c2(2)-c1(2)) / cfg.coarseDt;
    dExt   = min(1.20*Om_max/max(skySpd,1e-9), 0.35*(tEnd-tStart));

    T_prov = min((tEnd - tStart)/cfg.pathMul, ...
                 1.35 * max(T_est, cfg.horizonMul*max([t_dwell,t_slew,t_cent])));
    tStore = max(0, tStart - dExt);  dExt = tStart - tStore;
    span   = min(dExt + cfg.pathMul*T_prov, tEnd - tStore);
    tS     = 0 : cfg.storeDt : span;
    [rs, vs] = dacq_rk4J2(sr, sv, cfg.rk4Step, round(tStore/cfg.rk4Step));
    [AZs, ELs] = dacq_trackAzEl(rs, vs, cfg.epochJD + tStore/86400, tS, ...
                                cfg.station, 1.0);
    AZs = unwrapRows(AZs);
    Ks  = numel(tS);
    if verbose
        fprintf('  cloud spans %.1f s of centreline -> path extended by dExt = %.1f s\n', ...
                2*dExt/1.20, dExt);
    end

    % =====================================================================
    %  E. centreline / normal / width  (polynomial-smoothed)
    % =====================================================================
    gaz = median(AZs,1).';  gel = median(ELs,1).';
    sN  = 2*(tS(:) - tS(1))/(tS(end)-tS(1)) - 1;        % normalised time
    pAz = polyfit(sN, gaz, 6);   pEl = polyfit(sN, gel, 6);
    dAz = polyder(pAz);          dEl = polyder(pEl);
    ddAz= polyder(dAz);          ddEl= polyder(dEl);
    sc  = 2/(tS(end)-tS(1));                            % ds/dt

    GAMs  = [polyval(pAz,sN), polyval(pEl,sN)];
    tanAz = polyval(dAz,sN)*sc;  tanEl = polyval(dEl,sN)*sc;
    tx    = tanAz .* cos(GAMs(:,2));  ty = tanEl;
    nn    = hypot(tx,ty);  nn(nn==0) = 1;
    NHATs = [-ty./nn, tx./nn];

    % distribution width along the normal (2 sigma)
    Wds = zeros(Ks,1);
    for k = 1:Ks
        pn = (AZs(:,k)-GAMs(k,1)).*cos(GAMs(k,2))*NHATs(k,1) ...
           + (ELs(:,k)-GAMs(k,2))*NHATs(k,2);
        Wds(k) = 2*std(pn);
    end
    pW  = polyfit(sN, Wds, 4);  Wds = max(polyval(pW,sN), 1e-9);

    % =====================================================================
    %  F. dwell grid + Lagrange basis
    % =====================================================================
    pre = struct();
    pre.cfg=cfg; pre.pass=pass; pre.M=nPts; pre.dt=dt;
    pre.tStart=tStart; pre.tEnd=tEnd; pre.storeDt=cfg.storeDt;
    pre.dExt=dExt; pre.tStore=tStore; pre.skySpd=skySpd;
    pre.tS=tS; pre.AZs=AZs; pre.ELs=ELs; pre.Ks=Ks;
    pre.Om_max=Om_max; pre.nu=nu; pre.bottleneck=bottleneck;
    pre.N_spiral=N_spiral; pre.t_dwell=t_dwell; pre.t_slew=t_slew; pre.t_cent=t_cent;
    pre.b=b; pre.L_spiral=L_spiral;
    pre.polyAz=pAz; pre.polyEl=pEl; pre.polyW=pW; pre.sc=sc; pre.T_prov=T_prov;
    pre.GAMs=GAMs; pre.NHATs=NHATs; pre.Wds=Wds;

    N_prov  = floor(T_prov/dt) + 1;
    pre     = dacq_setgrid(pre, N_prov);

    % =====================================================================
    %  G. spiral baseline -> t_conv
    % =====================================================================
    [Psp, tTrav] = dacq_spiral(pre, cfg);
    res          = dacq_evaluate(Psp, pre, cfg);
    pre.spiral   = res; pre.Psp = Psp; pre.tTrav = tTrav;

    %  t_conv = the conventional spiral's FULL-TRAVERSE time (spec: "spiral
    %  완주 시간"). Reference coverage = what it has harvested by then. With
    %  c_pitch = 1.0 the footprints are tangent, so the spiral leaves
    %  interstitial gaps and never reaches 100 %; after the traverse it only
    %  creeps upward as the cloud drifts through a held pointing, which is
    %  not "scanning". Scoring every pattern by "time to reach the coverage
    %  the conventional method reaches, in the time it takes" is the fair
    %  comparison, and the same threshold is applied to all patterns.
    fr = min(cfg.completeFrac, floor(mean(res.t_hit <= tTrav)*100)/100);
    pre.frac_ref = fr;
    pre.t_conv   = tTrav;
    pre.t_max    = cfg.horizonMul * tTrav;
    N_max        = min(floor(pre.t_max/dt)+1, N_prov);
    pre          = dacq_setgrid(pre, N_max);

    % re-score the baseline on the final horizon so the baseline and every
    % optimised pattern are evaluated on exactly the same grid
    [pre.Psp, pre.tTrav] = dacq_spiral(pre, cfg);
    pre.spiral = dacq_evaluate(pre.Psp, pre, cfg);
    pre.spiral.t_reach = dacq_tfrac(pre.spiral, pre.frac_ref, pre.t_max);
    %  Like-for-like reference: the time the CONVENTIONAL scan needs to reach
    %  frac_ref. Comparing against its full-traverse time instead would flatter
    %  the optimiser, because most of the probability mass sits near the
    %  centre and the spiral collects it long before the traverse ends.
    pre.t_ref = pre.spiral.t_reach;
    if ~isfinite(pre.t_ref), pre.t_ref = pre.t_conv; end

    if verbose
        fprintf('  reference coverage = %.0f%% at the traverse (spiral %.1f%% by t_max)\n', ...
                100*pre.frac_ref, 100*pre.spiral.frac_scan);
        fprintf('  t_conv = %.2f s  t_max = %.2f s  N_max = %d  N_path = %d\n', ...
                pre.t_conv, pre.t_max, pre.N_max, pre.N_path);
    end
end

% =========================================================================
function A = unwrapRows(A)
% Remove 2*pi jumps along each row (azimuth continuity across the pass).
    for j = 1:size(A,1)
        A(j,:) = unwrap(A(j,:));
    end
    % keep all rows on the same branch as the first column median
    A = A - 2*pi*round((A(:,1) - median(A(:,1)))/(2*pi)) * ones(1,size(A,2));
end
