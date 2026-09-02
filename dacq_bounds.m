function [lo, hi] = dacq_bounds(pre, cfg)
%==========================================================================
% dacq_bounds  Box bounds of the 9 serpentine parameters (spec 4.2).
%   ln-amplitude nodes : 0.3*min(Wd) .. 3.0*max(Wd)
%   ln-accel   nodes   : 0.02*a_max  .. 1.00*a_max
%   drift multiplier   : 0.2 .. 3.0
%==========================================================================
    Wmin = min(pre.Wd);  Wmax = max(pre.Wd);
    lo = [repmat(log(0.3*Wmin),1,4), repmat(log(0.02*cfg.aMax),1,4), 0.2];
    hi = [repmat(log(3.0*Wmax),1,4), repmat(log(1.00*cfg.aMax),1,4), 3.0];
end
