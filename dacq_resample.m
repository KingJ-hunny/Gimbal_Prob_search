function pre2 = dacq_resample(pre, cfg, pass, nPts)
%==========================================================================
% dacq_resample  Same scan geometry, denser sampling of the same distribution.
%
%   Validating a tuned design by calling dacq_precompute again with more
%   points is WRONG: Om_max is a max over the sampled cloud, so a denser
%   sample enlarges it, which changes t_conv, the horizon, N_max, Wd and the
%   Lagrange basis. The design would then be scored on a different problem.
%
%   Here every geometric quantity is inherited from `pre` and only the sample
%   cloud is redrawn (same seed, so the smaller cloud is a subset), which is
%   the honest question: does a design tuned on M points still work when the
%   same distribution is sampled more densely?
%==========================================================================
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

    [rs, vs]   = dacq_rk4J2(sr, sv, cfg.rk4Step, round(pre.tStore/cfg.rk4Step));
    [AZs, ELs] = dacq_trackAzEl(rs, vs, cfg.epochJD + pre.tStore/86400, ...
                                pre.tS, cfg.station, 1.0);
    for j = 1:nPts, AZs(j,:) = unwrap(AZs(j,:)); end
    AZs = AZs - 2*pi*round((AZs(:,1)-median(AZs(:,1)))/(2*pi)) * ones(1,size(AZs,2));

    pre2 = pre;  pre2.AZs = AZs;  pre2.ELs = ELs;  pre2.M = nPts;
    pre2.spiral = dacq_evaluate(pre.Psp, pre2, cfg);
    pre2.t_ref  = dacq_tfrac(pre2.spiral, pre.frac_ref, pre.t_max);
    if ~isfinite(pre2.t_ref), pre2.t_ref = pre.t_conv; end
end
