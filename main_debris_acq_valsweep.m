%==========================================================================
% main_debris_acq_valsweep.m
%
%   Two corrections to the first experiment run:
%
%   (1) VALIDATION. Re-running dacq_precompute with 2000 points changes the
%       geometry (Om_max is a max over the sample, so it grows), which scores
%       the tuned design on a different problem. dacq_resample keeps the
%       geometry and only redraws the cloud.
%
%   (2) f_dwell SWEEP. It was run at (Q,R) = (1,1), where the optimum trades
%       coverage away to save control effort and often never reaches the 95 %
%       reference -- so the speedup was undefined. The sweep is a SPEED
%       question, so it belongs at (Q,R) = (1,0).
%
%   Both the baseline and the serpentine are scored with the same threshold.
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end
OUTDIR = 'figures';

cfg  = dacq_config();
pass = dacq_findpass(cfg);
pre  = dacq_precompute(cfg, pass, [], false);
load(fullfile(OUTDIR,'dacq_results.mat'));          % RES, SW, cfg, pass

% best speed-optimal (Q=1,R=0) design from the main run
[~, ib] = min(RES{2}.J);
xbest   = RES{2}.x(:,ib);
[~, best10] = dacq_objective(xbest, pre, cfg, 1, 0);
fprintf('\n=== reference (M = %d) ===\n', pre.M);
fprintf('  spiral     : t_ref %.2f s   J_x %.4f  J_u %.4f  frac %.3f\n', ...
        pre.t_ref, mean((pre.spiral.t_hit/pre.t_max).^2), ...
        mean((sqrt(sum(pre.spiral.u_hist.^2,2))/(cfg.aMax*pre.dt)).^2), ...
        pre.spiral.frac_scan);
fprintf('  serpentine : t_reach %s s  J_x %.4f  J_u %.4f  frac %.3f  speedup %s\n', ...
        fmtnum(best10.t_reach), best10.J_x, best10.J_u, best10.frac_scan, ...
        fmtnum(best10.speedup));

%% ---------------- (1) validation at 2000 points -----------------------
fprintf('\n=== validation: same geometry, %d points ===\n', cfg.nPointsVal);
preV = dacq_resample(pre, cfg, pass, cfg.nPointsVal);
resV = dacq_evaluate(dacq_pattern(best10.p, preV, cfg), preV, cfg);
tV   = dacq_tfrac(resV, pre.frac_ref, pre.t_max);
fprintf('  spiral     : frac %.3f  t_ref   %s s\n', preV.spiral.frac_scan, fmtnum(preV.t_ref));
fprintf('  serpentine : frac %.3f  t_reach %s s   speedup %s\n', ...
        resV.frac_scan, fmtnum(tV), fmtnum(preV.t_ref/tV));

fprintf('\n  coverage-time comparison on %d points\n', cfg.nPointsVal);
fprintf('  %8s | %10s | %12s | %8s\n','threshold','spiral [s]','serpentine [s]','ratio');
for fq = [0.5 0.8 0.9 0.95 0.99]
    a1 = dacq_tfrac(preV.spiral, fq, pre.t_max);
    a2 = dacq_tfrac(resV,        fq, pre.t_max);
    fprintf('  %7.0f%% | %10s | %12s | %8s\n', 100*fq, fmtnum(a1), fmtnum(a2), fmtnum(a1/a2));
end

%% ---------------- (2) f_dwell sweep at (Q,R) = (1,0) ------------------
fprintf('\n=== f_dwell sweep at (Q,R) = (1,0) ===\n');
FS = [100 250 500 1000 2000];  NT = 3;  BUD = 400;
SW = struct('f',[],'sp',[],'spq',[],'bn',{{}},'tconv',[],'tref',[],'nu',[],'jx',[],'jxsp',[]);
for i = 1:numel(FS)
    c2 = cfg; c2.fDwell = FS(i);
    p2 = dacq_precompute(c2, pass, cfg.nPoints, false);
    sp = []; jx = [];
    for tr = 1:NT
        X0 = dacq_initguess(p2, c2, 2, tr);
        f  = @(x) dacq_objective(x, p2, c2, 1, 0);
        xb = dacq_cmaes(f, X0(:,1), c2.sigma0, zeros(1,9), ones(1,9), BUD, 3000+tr, c2.popsize);
        [~, in] = dacq_objective(xb, p2, c2, 1, 0);
        if isfield(in,'speedup') && isfinite(in.speedup)
            sp(end+1) = in.speedup; jx(end+1) = in.J_x;
        end
    end
    if isempty(sp), sp = NaN; jx = NaN; end
    SW.f(end+1)=FS(i); SW.sp(end+1)=median(sp); SW.spq(:,end+1)=[qtl(sp,0.25);qtl(sp,0.75)];
    SW.bn{end+1}=p2.bottleneck; SW.tconv(end+1)=p2.t_conv; SW.tref(end+1)=p2.t_ref;
    SW.nu(end+1)=p2.nu; SW.jx(end+1)=median(jx);
    SW.jxsp(end+1)=mean((p2.spiral.t_hit/p2.t_max).^2);
    fprintf('  f=%5d Hz | %-11s | t_conv %6.2f | t_ref %6.2f | J_x spiral %.4f vs serp %.4f | speedup %.2f\n', ...
            FS(i), p2.bottleneck, p2.t_conv, p2.t_ref, SW.jxsp(end), SW.jx(end), SW.sp(end));
end
dacq_plots('sweep', pre, cfg, OUTDIR, SW);
dacq_plots('cumulative', preV, cfg, OUTDIR, resV);
dacq_plots('control',    pre,  cfg, OUTDIR, best10);
dacq_plots('pattern',    pre,  cfg, OUTDIR, best10);
dacq_plots('profiles',   pre,  cfg, OUTDIR, best10);

save('-mat7-binary', fullfile(OUTDIR,'dacq_valsweep.mat'), 'SW','best10','tV');
fprintf('\nDone.\n');
