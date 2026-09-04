%==========================================================================
% main_debris_acq_skytime.m
%   Dwell points of the spiral and the optimised serpentine on the absolute
%   az/el sky, drawn together with the mean-orbit profile and coloured by
%   time. Uses the (Q,R) = (1,0) design from main_debris_acq_refair.
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end
OUTDIR='figures';
cfg=dacq_config(); pass=dacq_findpass(cfg); pre=dacq_precompute(cfg,pass,[],false);
load(fullfile(OUTDIR,'dacq_refair.mat'));            % OUT
[~,ib]=min(OUT{1}.J);
[~,best]=dacq_objective(OUT{1}.x(:,ib), pre, cfg, 1, 0);
dacq_plots('skytime', pre, cfg, OUTDIR, best);
fprintf('mean orbit over the scan: az %.2f -> %.2f deg, el %.2f -> %.2f deg in %.2f s\n', ...
   rad2deg(pre.GAM(1,1)), rad2deg(pre.GAM(end,1)), ...
   rad2deg(pre.GAM(1,2)), rad2deg(pre.GAM(end,2)), pre.t_max);
fprintf('scan offset extent: +/- %.0f arcsec (Om_max)\n', pre.Om_max*206264.8);
fprintf('Done.\n');
