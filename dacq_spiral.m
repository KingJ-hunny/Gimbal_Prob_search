function [P, tTrav] = dacq_spiral(pre, cfg)
%==========================================================================
% dacq_spiral  Conventional Archimedean baseline, ACTUATOR-FEASIBLE.
%
%   This wraps dacq_spiral2 with the policy that makes the baseline a fair
%   opponent for the optimised serpentine:
%
%     profile 'proper'  the speed profile bounds the TOTAL acceleration
%                       (tangential and centripetal together) and is marched
%                       in time on the dwell grid, so the SAMPLED path is
%                       feasible. The earlier version capped only the
%                       centripetal term via v <= sqrt(a*R); that leaves the
%                       tangential term free, and the traverse ended by
%                       setting the speed to zero in one sample -- a velocity
%                       step of 145-185 deg/s^2 against a 5 deg/s^2 limit.
%                       The serpentine was always held to the limit by the
%                       objective's feasibility test, but the spiral never
%                       went through that test, so the two were not being
%                       held to the same standard.
%     accPolicy 'inst'  the same acceleration budget the serpentine gets
%                       (a_avail), de-rated by the instantaneous cos(el) that
%                       converts an on-sky offset into azimuth-axis motion.
%     useDwell true     the dwell-spacing cap 2*r_beam*f_dwell is kept: a
%                       conventional spiral is designed to tile without gaps.
%
%   The previous behaviour is preserved in dacq_spiral_naive.m for comparison.
%==========================================================================
    opt = struct('accPolicy','inst', 'velFac',1.00, ...
                 'useDwell',true,    'profile','proper');
    [P, tTrav] = dacq_spiral2(pre, cfg, opt);
end
