%==========================================================================
% main_debris_acq_figs.m
%   Regenerate the Stage-7 diagnostics. The f_dwell sweep numbers below are
%   the medians measured by main_debris_acq_valsweep (3 trials, 400 evals per
%   point, (Q,R) = (1,0)); that run crashed in its final plotting call before
%   it could save, so they are transcribed here rather than re-optimised.
%   Per-trial spread was not retained, so no IQR band is drawn.
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end
OUTDIR = 'figures';

cfg  = dacq_config();
pass = dacq_findpass(cfg);
pre  = dacq_precompute(cfg, pass, [], false);
load(fullfile(OUTDIR,'dacq_results.mat'));

[~, ib] = min(RES{2}.J);
[~, best10] = dacq_objective(RES{2}.x(:,ib), pre, cfg, 1, 0);
preV = dacq_resample(pre, cfg, pass, cfg.nPointsVal);
resV = dacq_evaluate(dacq_pattern(best10.p, preV, cfg), preV, cfg);
resV.t_reach = dacq_tfrac(resV, pre.frac_ref, pre.t_max);

SW.f     = [100 250 500 1000 2000];
SW.sp    = [0.94 0.89 0.88 0.90 0.86];
SW.spq   = nan(2,5);
SW.tconv = [16.36 16.36 16.36 16.36 16.36];
SW.tref  = [9.79 9.20 9.20 9.20 9.20];
SW.jxsp  = [0.0353 0.0274 0.0273 0.0272 0.0272];
SW.jx    = [0.1016 0.1101 0.1070 0.1056 0.1063];
SW.nu    = repmat(pre.nu,1,5);
SW.bn    = {'centripetal','centripetal','centripetal','centripetal','centripetal'};

dacq_plots('centerline', pre,  cfg, OUTDIR);
dacq_plots('spiral',     pre,  cfg, OUTDIR);
dacq_plots('pattern',    pre,  cfg, OUTDIR, best10);
dacq_plots('profiles',   pre,  cfg, OUTDIR, best10);
dacq_plots('thist',      pre,  cfg, OUTDIR, best10);
dacq_plots('control',    pre,  cfg, OUTDIR, best10);
dacq_plots('cumulative', preV, cfg, OUTDIR, resV);
dacq_plots('sweep',      pre,  cfg, OUTDIR, SW);
save('-mat7-binary', fullfile(OUTDIR,'dacq_valsweep.mat'), 'SW','best10','resV');
fprintf('serpentine (1,0): frac %.3f  t_reach %.2f  speedup %.3f  J_x %.4f  J_u %.4f\n', ...
        best10.frac_scan, best10.t_reach, best10.speedup, best10.J_x, best10.J_u);
fprintf('validation 2000 pts: spiral %.2f s, serpentine %.2f s, speedup %.2f\n', ...
        preV.t_ref, resV.t_reach, preV.t_ref/resV.t_reach);
fprintf('Done.\n');
