function res = dacq_evaluate(P, pre, cfg)
%==========================================================================
% dacq_evaluate  Binary-detection scoring of a dwell sequence.
%
%   Sample point j is detected the first time a dwell lands within r_beam of
%   it; detected points are harvested (removed). Because each point's first
%   hit is independent of the others, the harvest is an efficiency device
%   only -- it does not change t_hit -- so the loop runs over CHUNKS of
%   epochs with the surviving points fully vectorised inside each chunk.
%
%   Given:  P    N x 2 dwell directions [az el] (rad)
%           pre  precomputed layer-0 struct;  cfg  config
%   Return: res  .t_hit .frac_scan .t_complete .tSorted .u_hist .N_used
%==========================================================================
    M  = pre.M;  N = min(size(P,1), pre.N_max);  dt = pre.dt;
    tmax = pre.t_max;  r2 = cfg.rBeam^2;
    sdt = pre.storeDt;  Ks = pre.Ks;

    t_hit  = tmax * ones(M,1);
    active = true(M,1);
    tD     = pre.tD;
    B      = 1024;
    N_used = N;  reason = 'horizon';

    k0 = 1;
    while k0 <= N
        if ~any(active)
            N_used = k0 - 1;  reason = 'complete';  break;
        end
        k1  = min(k0+B-1, N);
        idx = find(active);

        % ---- interpolate the cloud onto this chunk of dwell epochs -----
        %  the storage grid starts dExt before t_start (see dacq_setgrid)
        tt = tD(k0:k1).' + pre.dExt;                       % 1 x Bc
        ii = min(max(floor(tt/sdt)+1, 1), Ks-1);
        ff = (tt - (ii-1)*sdt) / sdt;
        A1 = pre.AZs(idx,ii);  A2 = pre.AZs(idx,ii+1);
        E1 = pre.ELs(idx,ii);  E2 = pre.ELs(idx,ii+1);
        AZb = A1 + (A2-A1).*ff;
        ELb = E1 + (E2-E1).*ff;

        % ---- squared angular distance (cos(el) metric, no arccos) ------
        daz = P(k0:k1,1).' - AZb;
        del = P(k0:k1,2).' - ELb;
        d2  = (daz .* cos(ELb)).^2 + del.^2;

        % ---- first hit per surviving point ----------------------------
        [mx, col] = max(d2 <= r2, [], 2);       % max finds the FIRST 1
        h = mx > 0;
        if any(h)
            rows = idx(h);
            t_hit(rows)  = tD(k0 + col(h) - 1);
            active(rows) = false;
        end
        k0 = k1 + 1;
    end

    % ---- control effort: u_k = w(k-1->k) - w(k-2->k-1) -----------------
    nu_ = max(N_used-2, 0);
    if nu_ > 0
        U = diff(P(1:min(N_used,size(P,1)),:), 2, 1) / dt;
        res.u_hist = U(1:min(nu_,size(U,1)),:);
    else
        res.u_hist = zeros(0,2);
    end

    res.t_hit      = t_hit;
    res.tSorted    = sort(t_hit);
    res.frac_scan  = mean(t_hit < tmax);
    res.N_used     = N_used;
    res.reason     = reason;
    if strcmp(reason,'complete'), res.t_complete = N_used*dt; else, res.t_complete = tmax; end
end
