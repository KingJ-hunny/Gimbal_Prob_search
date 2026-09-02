function pre = dacq_setgrid(pre, N)
%==========================================================================
% dacq_setgrid  Build the PATH grid and the EVALUATION grid for N dwells.
%
%   Why two grids. A sample point displaced along-track by s is seen, at any
%   instant, in the direction the median occupies s/v_orbit seconds later, so
%   the cloud lies ALONG the centreline path, spanning +/- dExt seconds of it.
%   The spec drives the pointing with idx = v*k starting at path index 0, a
%   monotone sweep; if the path began at t_start it could only ever reach the
%   LEADING half of the cloud (~50 % ceiling, no matter what v is).
%   So the path is extended by dExt at BOTH ends: path index 0 is where the
%   trailing edge sits at t_start, and a single monotone sweep with v > 1 can
%   cover the whole cloud. The 9-parameter structure is untouched.
%
%     path grid  : N_path = N + 2*nExt dwell steps, used only as GEOMETRY
%                  for the pointing to traverse (GAMp / NHATp).
%     eval grid  : N steps of real time from t_start (tD), carrying Wd,
%                  a_avail and the Lagrange basis, plus the centreline as
%                  the cloud actually moves (GAMe / NHATe) for the baseline.
%==========================================================================
    cfg = pre.cfg;  dt = pre.dt;  T = pre.tS(end);
    nExt  = round(pre.dExt/dt);
    % long enough that idx = v*(N-1) never saturates for v up to pathMul
    Npath = min(floor(T/dt) + 1, ceil(cfg.pathMul*(N-1)) + 1);
    Npath = max(Npath, N + nExt);

    dAz = polyder(pre.polyAz);  ddAz = polyder(dAz);
    dEl = polyder(pre.polyEl);  ddEl = polyder(dEl);
    sc  = 2/T;

    % ---- path grid (geometry the pointing may traverse) ----------------
    tP   = (0:Npath-1).' * dt;
    sP   = 2*tP/T - 1;
    GAMp = [polyval(pre.polyAz,sP), polyval(pre.polyEl,sP)];
    txP  = polyval(dAz,sP)*sc .* cos(GAMp(:,2));
    tyP  = polyval(dEl,sP)*sc;
    nP   = hypot(txP,tyP);  nP(nP==0) = 1;
    NHATp= [-tyP./nP, txP./nP];

    % ---- evaluation grid (real time, starting at t_start) --------------
    ke   = nExt + (1:N).';
    tE   = tP(ke);
    sE   = 2*tE/T - 1;
    GAMe = GAMp(ke,:);   NHATe = NHATp(ke,:);

    Wd   = max(polyval(pre.polyW, sE), 1e-9);
    ffAz = polyval(ddAz,sE)*sc^2;
    ffEl = polyval(ddEl,sE)*sc^2;
    aAvail = cfg.rho * max(cfg.aMax - max(abs(ffAz),abs(ffEl)), 0.05*cfg.aMax);

    nodes = [-1.0 -0.5 0.5 1.0];
    sN    = 2*(0:N-1).'/(N-1) - 1;              % pattern's own normalised step
    LEG   = ones(N,4);
    for i = 1:4
        Li = ones(N,1);
        for j = 1:4
            if j ~= i, Li = Li .* (sN - nodes(j)) / (nodes(i) - nodes(j)); end
        end
        LEG(:,i) = Li;
    end

    pre.N_max   = N;      pre.N_path = Npath;   pre.nExt = nExt;
    pre.tD      = (0:N-1).' * dt;               % time since t_start
    pre.t_max   = (N-1)*dt;
    pre.GAMp    = GAMp;   pre.NHATp = NHATp;
    pre.GAM     = GAMe;   pre.NHAT  = NHATe;    % centreline in real time
    pre.Wd      = Wd;     pre.a_avail = aAvail; pre.LEG = LEG;
    pre.nodes   = nodes;
end
