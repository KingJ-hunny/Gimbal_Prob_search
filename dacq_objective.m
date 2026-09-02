function [J, info] = dacq_objective(x, pre, cfg, Q, R)
%==========================================================================
% dacq_objective  J = Q*J_x + R*J_u for a normalised parameter vector.
%
%   Feasibility is a hard constraint, but the returned value CARRIES THE
%   VIOLATION MAGNITUDE: returning a flat +1e6 would make every infeasible
%   candidate tie, and CMA-ES would lose the direction back to the feasible
%   set.
%
%   Both cost terms are normalised to [0,1] before weighting -- t_hit is in
%   seconds (tens) and u is in rad/s (1e-4), so an unnormalised (1,1) would
%   be indistinguishable from (1,0).
%==========================================================================
    [lo, hi] = dacq_bounds(pre, cfg);
    p  = lo + x(:).' .* (hi - lo);
    dt = pre.dt;

    [P, A, a, OM] = dacq_pattern(p, pre, cfg);

    vel  = diff(P,1,1)/dt;
    acc  = diff(P,2,1)/dt^2;
    viol = max(0, max(abs(vel(:))) - cfg.wMax)/cfg.wMax ...
         + max(0, max(abs(acc(:))) - cfg.aMax)/cfg.aMax ...
         + max(0, deg2rad(cfg.elMaskDeg) - min(P(:,2)))/deg2rad(cfg.elMaskDeg);
    if viol > 0
        J = 1e6 + viol;
        if nargout > 1, info = struct('viol',viol,'feasible',false); end
        return;
    end

    res = dacq_evaluate(P, pre, cfg);

    J_x = mean( (res.t_hit / pre.t_max).^2 );
    if isempty(res.u_hist)
        J_u = 0;
    else
        J_u = mean( (sqrt(sum(res.u_hist.^2,2)) / (cfg.aMax*dt)).^2 );
    end
    J = Q*J_x + R*J_u;

    if nargout > 1
        info          = res;
        info.J_x      = J_x;
        info.J_u      = J_u;
        info.J        = J;
        info.viol     = 0;
        info.feasible = true;
        info.p        = p;
        info.P        = P;
        info.A        = A;
        info.a        = a;
        info.OM       = OM;
        info.t_reach  = dacq_tfrac(res, pre.frac_ref, pre.t_max);
        info.speedup  = pre.t_ref / info.t_reach;
    end
end
