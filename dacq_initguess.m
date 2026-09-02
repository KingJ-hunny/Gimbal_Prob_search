function X0 = dacq_initguess(pre, cfg, nRandom, seed)
%==========================================================================
% dacq_initguess  Physically distinct starting basins (spec 6.1), returned
%                 as normalised [0,1] column vectors.
%
%   (1) constant acceleration : zeta flat  -- actuator used evenly
%   (2) constant speed        : zeta = 2 ln c - eta  -- Archimedean-equivalent
%   (3) constant frequency    : zeta = 2 ln Om0 + eta
%   (4) uniform               : eta and zeta both flat
%   plus nRandom uniform random starts.
%==========================================================================
    if nargin < 3, nRandom = 2; end
    if nargin < 4, seed = 0; end
    [lo, hi] = dacq_bounds(pre, cfg);
    N = pre.N_max;

    kn  = max(1, min(N, round(((pre.nodes+1)/2)*(N-1)) + 1));
    Wn  = pre.Wd(kn).';                       % width at the 4 nodes
    eta = log(1.2*Wn);
    Am  = mean(1.2*Wn);
    a0  = 0.5*cfg.aMax;
    c0  = sqrt(a0*Am);                        % transverse speed scale
    Om0 = sqrt(a0/Am);                        % frequency scale

    % Drift multiplier. v = 1 (the spec's value) just rides the trailing edge
    % of the cloud and scans nothing, so the physical starts use the GAP-FREE
    % maximum instead: the pointing may advance along the path by at most one
    % beam diameter per half oscillation, (v-1)*skySpd*(pi/Omega) <= 2*r_beam.
    Om0g = sqrt(a0/Am);
    vgap = 1 + 2*cfg.rBeam*Om0g / (pi * max(pre.skySpd,1e-9));
    vs   = min(max(vgap, 0.2), 3.0);

    g1 = [eta,                    repmat(log(a0),1,4),        vs ];
    g2 = [eta,                    2*log(c0) - eta,            vs ];
    g3 = [eta,                    2*log(Om0) + eta,           vs ];
    g4 = [repmat(log(mean(pre.Wd)),1,4), repmat(log(a0),1,4), 1.0];

    G = [g1; g2; g3; g4];
    X0 = zeros(9, size(G,1) + nRandom);
    for i = 1:size(G,1)
        X0(:,i) = min(max((G(i,:) - lo) ./ (hi - lo), 0), 1).';
    end
    rand('seed', seed + 777);
    for i = 1:nRandom
        X0(:, size(G,1)+i) = rand(9,1);
    end
end
