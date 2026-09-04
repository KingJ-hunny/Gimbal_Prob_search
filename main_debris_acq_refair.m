%==========================================================================
% main_debris_acq_refair.m
%   Re-run the (Q,R) comparison against the ACTUATOR-FEASIBLE baseline.
%   The original dacq_spiral capped only the centripetal acceleration and
%   stopped dead at the end of its traverse (145-185 deg/s^2 against a
%   5 deg/s^2 limit), and unlike the serpentine it never passed through the
%   objective's feasibility test. Both patterns are now held to the same
%   actuator limits.
%==========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), try, graphics_toolkit('gnuplot'); catch, end; end
OUTDIR='figures';
cfg=dacq_config(); pass=dacq_findpass(cfg); pre=dacq_precompute(cfg,pass);
dt=pre.dt;
vv=diff(pre.Psp,1,1)/dt; aa=diff(pre.Psp,2,1)/dt^2;
spJx=mean((pre.spiral.t_hit/pre.t_max).^2);
spJu=mean((sqrt(sum(pre.spiral.u_hist.^2,2))/(cfg.aMax*dt)).^2);
fprintf('\nBASELINE (feasible): max|v| %.2f deg/s  max|a| %.2f deg/s^2  J_x %.4f  J_u %.4f  t_ref %.2f s\n',...
   rad2deg(max(abs(vv(:)))), rad2deg(max(abs(aa(:)))), spJx, spJu, pre.t_ref);

QR=[1 0; 1 1]; nm={'Q=1,R=0','Q=1,R=1'}; OUT=cell(1,2);
for q=1:2
  Q=QR(q,1); R=QR(q,2); T=struct('J',[],'frac',[],'tr',[],'sp',[],'jx',[],'ju',[],'x',[]);
  fprintf('\n=== (%s) 4 trials x 2 starts x 700 evals ===\n', nm{q});
  for tr=1:4
    X0=dacq_initguess(pre,cfg,4,tr); pick=[1, 4+mod(tr,4)+1];
    bJ=inf; bx=[];
    for s=pick
      f=@(x) dacq_objective(x,pre,cfg,Q,R);
      [xb,fb]=dacq_cmaes(f,X0(:,s),cfg.sigma0,zeros(1,9),ones(1,9),700,1000*tr+s,cfg.popsize);
      if fb<bJ, bJ=fb; bx=xb; end
    end
    [~,in]=dacq_objective(bx,pre,cfg,Q,R);
    T.J(end+1)=bJ; T.x(:,end+1)=bx;
    T.frac(end+1)=in.frac_scan; T.tr(end+1)=in.t_reach; T.sp(end+1)=in.speedup;
    T.jx(end+1)=in.J_x; T.ju(end+1)=in.J_u;
    fprintf('  trial %d: J=%.4f frac=%.3f t_reach=%s speedup=%s\n',tr,bJ,in.frac_scan,fmtnum(in.t_reach),fmtnum(in.speedup));
  end
  OUT{q}=T;
  fprintf('  --- median [IQR] ---\n');
  repIQR('frac_scan',T.frac); repIQR('t_reach [s]',T.tr); repIQR('speedup',T.sp);
  repIQR('J_x',T.jx); repIQR('J_u',T.ju);
end
[~,ib]=min(OUT{1}.J); [~,b10]=dacq_objective(OUT{1}.x(:,ib),pre,cfg,1,0);
dacq_plots('cumulative',pre,cfg,OUTDIR,b10);
dacq_plots('control',   pre,cfg,OUTDIR,b10);
dacq_plots('thist',     pre,cfg,OUTDIR,b10);
dacq_plots('spiral',    pre,cfg,OUTDIR);
save('-mat7-binary',fullfile(OUTDIR,'dacq_refair.mat'),'OUT','spJx','spJu');
fprintf('\nDone.\n');
