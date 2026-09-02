function [P, A, a, OM, PHI] = dacq_pattern(p, pre, cfg)
%==========================================================================
% dacq_pattern  Serpentine dwell sequence from the 9 design parameters.
%
%   p = [eta(1:4) zta(1:4) v]
%     eta : ln-amplitude at the 4 Lagrange nodes   -> A(s) = exp(LEG*eta)
%     zta : ln-acceleration budget at those nodes  -> a(s) = exp(LEG*zta)
%     v   : centreline drift-speed multiplier
%
%   Working in the log domain makes A and a positive by construction, and
%   binding Omega through A*Omega^2 = a turns the acceleration limit into a
%   plain box constraint on a. The phase is obtained by TRAPEZOIDAL
%   INTEGRATION of Omega, not as Omega*t, so the pattern stays continuous
%   however fast Omega varies.
%==========================================================================
    N = pre.N_max;  dt = pre.dt;
    eta = p(1:4); eta = eta(:);
    zta = p(5:8); zta = zta(:);
    v   = p(9);

    A  = exp(pre.LEG * eta);
    a  = exp(pre.LEG * zta);
    a  = min(a, pre.a_avail);            % never spend more than is available

    OM  = sqrt(a ./ A) * dt;             % [rad per step]
    PHI = [0; cumsum(0.5*(OM(1:end-1) + OM(2:end)))];

    % ---- centreline PATH sampled at the drifting index -----------------
    %  Index 0 of the path is where the trailing edge of the cloud sits at
    %  t_start, so a monotone sweep with v > 1 crosses the whole cloud.
    Np  = pre.N_path;
    idx = min(max(v*(0:N-1).', 0), Np-1);
    i0  = floor(idx);  fr = idx - i0;  i0 = i0 + 1;
    i1  = min(i0+1, Np);
    GAMk  = pre.GAMp(i0,:)  + (pre.GAMp(i1,:)  - pre.GAMp(i0,:) ).*fr;
    NHATk = pre.NHATp(i0,:) + (pre.NHATp(i1,:) - pre.NHATp(i0,:)).*fr;

    % ---- dwell direction: tangent-plane offset mapped back to az/el ----
    n    = A .* sin(PHI);
    P_el = GAMk(:,2) + n .* NHATk(:,2);
    P_az = GAMk(:,1) + n .* NHATk(:,1) ./ cos(P_el);
    P    = [P_az, P_el];
end
