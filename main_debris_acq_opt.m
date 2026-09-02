%==========================================================================
% main_debris_acq_opt.m
%
%   Debris laser-ranging acquisition: CMA-ES optimisation of a 9-parameter
%   serpentine scan against a conventional Archimedean-spiral baseline.
%
%   Pipeline (implementation spec sections 3-7):
%     LAYER 0   dacq_findpass -> dacq_precompute   (cloud, centreline, t_conv)
%     baseline  dacq_spiral   + dacq_evaluate
%     design    dacq_pattern  (9 params) -> dacq_objective (binary detection)
%     solver    dacq_cmaes    (Hansen CMA-ES, implemented here: Octave has no
%                              `cma` package, and the binary detection test
%                              makes the objective piecewise constant, so a
%                              gradient method is not available anyway)
%
%   Reported metric:  speedup = t_ref / t_reach, both being the time to harvest
%   the SAME fraction (frac_ref) of the probability cloud -- the conventional
%   spiral's time in the numerator. Comparing against the spiral's full
%   traverse time instead would flatter the optimiser, because most of the
%   probability mass sits near the centreline and the centre-out spiral
%   collects it long before its traverse ends.
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end
OUTDIR = 'figures'; if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

% ---- experiment size (see cfg SCALE block for why this is not the spec's
%      3 x 30 trials x 8 starts x 5000 evals matrix) ----------------------
BUDGET   = 700;      % CMA-ES evaluations per start
NTRIAL   = 4;        % repeats per (Q,R)
NSTART   = 2;        % starts per trial (1 physical + 1 random)
SWEEP_B  = 400;      % budget for the f_dwell sweep
SWEEP_T  = 3;        % trials for the f_dwell sweep
FSWEEP   = [100 250 500 1000 2000];

cfg  = dacq_config();
fprintf('\n================ LAYER 0 ================\n');
pass = dacq_findpass(cfg);
fprintf('  pass: peak el %.1f deg at %.2f h after epoch (RAAN %.0f, M0 %.0f)\n', ...
        pass.elPeakDeg, pass.tPeak/3600, pass.raanDeg, pass.M0Deg);
pre  = dacq_precompute(cfg, pass);

fprintf('\n  GATE nu           : %.1f  (%s)\n', pre.nu, ...
        ternary(pre.nu>5 && pre.nu<500,'PASS','FAIL'));
fprintf('  GATE bottleneck   : %s\n', pre.bottleneck);
fprintf('  GATE t_conv       : %.2f s  (%s)\n', pre.t_conv, ...
        ternary(pre.t_conv>1 && pre.t_conv<500,'plausible','CHECK'));
fprintf('  aspect ratio      : Om_max/mean(Wd) = %.1f  (Wd %.0f..%.0f arcsec)\n', ...
        pre.Om_max/mean(pre.Wd), min(pre.Wd)*206264.8, max(pre.Wd)*206264.8);
fprintf('  baseline spiral   : frac %.3f by t_max, reference coverage %.0f%%\n', ...
        pre.spiral.frac_scan, 100*pre.frac_ref);

%  Score the BASELINE under the very objective being optimised. If the
%  spiral's J_x is already lower than anything the serpentine family reaches,
%  the honest conclusion is that the conventional scan wins on speed.
spJx = mean((pre.spiral.t_hit/pre.t_max).^2);
spJu = mean((sqrt(sum(pre.spiral.u_hist.^2,2))/(cfg.aMax*pre.dt)).^2);
fprintf('  baseline scored as a candidate: J_x = %.4f, J_u = %.4f\n', spJx, spJu);
fprintf('  t_ref (spiral to %.0f%%) = %.2f s ; full traverse t_conv = %.2f s\n', ...
        100*pre.frac_ref, pre.t_ref, pre.t_conv);

dacq_plots('centerline', pre, cfg, OUTDIR);
dacq_plots('spiral',     pre, cfg, OUTDIR);

%% ===================== test matrix over (Q,R) ==========================
QR   = [1 1; 1 0; 0 1];
names= {'Q=1,R=1','Q=1,R=0','Q=0,R=1'};
RES  = cell(1,3);
for q = 1:3
    Q = QR(q,1); R = QR(q,2);
    nt = NTRIAL; nb = BUDGET; ns = NSTART;
    if Q == 0, nt = 1; nb = 150; ns = 1; end        % degenerate sanity check
    fprintf('\n============ (%s)  %d trials x %d starts x %d evals ============\n', ...
            names{q}, nt, ns, nb);
    T = struct('J',[],'frac',[],'treach',[],'speedup',[],'Jx',[],'Ju',[],'x',[]);
    for tr = 1:nt
        X0 = dacq_initguess(pre, cfg, 4, tr);
        pick = [1, 4+mod(tr,4)+1];  pick = pick(1:ns);
        bJ = inf; bx = [];
        for s = pick
            f = @(x) dacq_objective(x, pre, cfg, Q, R);
            [xb, fb] = dacq_cmaes(f, X0(:,s), cfg.sigma0, zeros(1,9), ones(1,9), ...
                                  nb, 1000*tr+s, cfg.popsize);
            if fb < bJ, bJ = fb; bx = xb; end
        end
        [~, in] = dacq_objective(bx, pre, cfg, Q, R);
        T.J(end+1)=bJ; T.x(:,end+1)=bx;
        if isfield(in,'frac_scan')
            T.frac(end+1)=in.frac_scan; T.treach(end+1)=in.t_reach;
            T.speedup(end+1)=in.speedup; T.Jx(end+1)=in.J_x; T.Ju(end+1)=in.J_u;
        else
            T.frac(end+1)=0; T.treach(end+1)=NaN; T.speedup(end+1)=NaN;
            T.Jx(end+1)=NaN; T.Ju(end+1)=NaN;
        end
        fprintf('   trial %d: J=%.4f frac=%.3f t_reach=%s speedup=%s\n', tr, bJ, ...
                T.frac(end), fmtnum(T.treach(end)), fmtnum(T.speedup(end)));
    end
    RES{q} = T;
    fprintf('   --- median [IQR] over %d trials ---\n', nt);
    repIQR('frac_scan', T.frac); repIQR('t_reach [s]', T.treach);
    repIQR('speedup',   T.speedup); repIQR('J_x', T.Jx); repIQR('J_u', T.Ju);
end

%% ---- diagnostics for the main operating point (1,1) -------------------
[~, ib] = min(RES{1}.J);
[~, best11] = dacq_objective(RES{1}.x(:,ib), pre, cfg, 1, 1);
dacq_plots('pattern',    pre, cfg, OUTDIR, best11);
dacq_plots('profiles',   pre, cfg, OUTDIR, best11);
dacq_plots('thist',      pre, cfg, OUTDIR, best11);
dacq_plots('cumulative', pre, cfg, OUTDIR, best11);
dacq_plots('control',    pre, cfg, OUTDIR, best11);

%% ---- validation on the full 2000-point cloud --------------------------
fprintf('\n============ validation on %d points ============\n', cfg.nPointsVal);
preV = dacq_precompute(cfg, pass, cfg.nPointsVal, false);
resV = dacq_evaluate(dacq_pattern(best11.p, preV, cfg), preV, cfg);
tV   = dacq_tfrac(resV, preV.frac_ref, preV.t_max);
fprintf('  spiral    : t_ref %.2f s (traverse %.2f s), frac %.3f\n', ...
        preV.t_ref, preV.t_conv, preV.spiral.frac_scan);
fprintf('  serpentine: frac %.3f, t_reach %s s, speedup %s\n', ...
        resV.frac_scan, fmtnum(tV), fmtnum(preV.t_ref/tV));

%% ===================== f_dwell sweep ==================================
fprintf('\n============ f_dwell sweep ============\n');
SW = struct('f',[],'sp',[],'spq',[],'bn',{{}},'tconv',[],'nu',[]);
for i = 1:numel(FSWEEP)
    c2 = cfg; c2.fDwell = FSWEEP(i);
    p2 = dacq_precompute(c2, pass, cfg.nPoints, false);
    sp = [];
    for tr = 1:SWEEP_T
        X0 = dacq_initguess(p2, c2, 2, tr);
        f  = @(x) dacq_objective(x, p2, c2, 1, 1);
        [xb, ~] = dacq_cmaes(f, X0(:,1), c2.sigma0, zeros(1,9), ones(1,9), ...
                             SWEEP_B, 2000+tr, c2.popsize);
        [~, in] = dacq_objective(xb, p2, c2, 1, 1);
        if isfield(in,'speedup') && isfinite(in.speedup), sp(end+1) = in.speedup; end
    end
    if isempty(sp), sp = NaN; end
    SW.f(end+1)=FSWEEP(i); SW.sp(end+1)=median(sp);
    SW.spq(:,end+1)=[qtl(sp,0.25); qtl(sp,0.75)];
    SW.bn{end+1}=p2.bottleneck; SW.tconv(end+1)=p2.t_conv; SW.nu(end+1)=p2.nu;
    fprintf('  f=%5d Hz | bottleneck %-11s | t_conv %6.2f s | speedup median %.2f\n', ...
            FSWEEP(i), p2.bottleneck, p2.t_conv, SW.sp(end));
end
dacq_plots('sweep', pre, cfg, OUTDIR, SW);

save('-mat7-binary', fullfile(OUTDIR,'dacq_results.mat'), 'RES','SW','cfg','pass');
fprintf('\nFigures + dacq_results.mat written to %s/\n', OUTDIR);
fprintf('Done.\n');
