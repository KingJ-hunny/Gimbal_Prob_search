function [xbest, fbest, out] = dacq_cmaes(fun, x0, sigma0, lb, ub, maxfe, seed, lambda)
%==========================================================================
% dacq_cmaes  (mu/mu_w, lambda)-CMA-ES with box handling by clipping.
%
%   Hansen's standard formulation. Octave has no `cma` package, so this is a
%   direct implementation. Gradient-free is mandatory here: the binary
%   detection test makes the objective piecewise constant.
%
%   sigma0 must stay large (0.3 of the normalised box). If it is small every
%   member of a generation lands in the same constant piece, all ranks tie,
%   and the algorithm gets no information.
%
%   Given:  fun  handle f(x) with x a normalised 9-vector
%           x0   start (column), sigma0 step, lb/ub bounds (row)
%           maxfe evaluation budget, seed, lambda population size
%   Return: xbest, fbest, out (.fhist best-so-far per generation, .nfe)
%==========================================================================
    n  = numel(x0);
    x0 = x0(:);  lb = lb(:);  ub = ub(:);
    if nargin < 8 || isempty(lambda), lambda = 4 + floor(3*log(n)); end
    mu = floor(lambda/2);
    w  = log(mu+0.5) - log(1:mu).';  w = w/sum(w);
    mueff = 1/sum(w.^2);

    cc    = (4 + mueff/n) / (n + 4 + 2*mueff/n);
    cs    = (mueff + 2) / (n + mueff + 5);
    c1    = 2 / ((n+1.3)^2 + mueff);
    cmu   = min(1-c1, 2*(mueff - 2 + 1/mueff) / ((n+2)^2 + mueff));
    damps = 1 + 2*max(0, sqrt((mueff-1)/(n+1)) - 1) + cs;
    chiN  = sqrt(n)*(1 - 1/(4*n) + 1/(21*n^2));

    randn('seed', seed);  rand('seed', seed);
    xmean = x0;  sigma = sigma0;
    pc = zeros(n,1);  ps = zeros(n,1);
    B  = eye(n);  D = ones(n,1);  C = eye(n);  invsqrtC = eye(n);

    % score the start itself: CMA-ES only ever evaluates SAMPLED points, so
    % without this a good x0 can be lost entirely when sigma0 is large.
    xbest = min(max(x0,lb),ub);  fbest = fun(xbest);  nfe = 1;
    fhist = [];  gen = 0;
    while nfe < maxfe
        gen = gen + 1;
        arx = zeros(n, lambda);  arf = zeros(1, lambda);
        for k = 1:lambda
            y  = B * (D .* randn(n,1));
            xk = min(max(xmean + sigma*y, lb), ub);     % clip into the box
            arx(:,k) = xk;
            arf(k)   = fun(xk);
            nfe = nfe + 1;
        end
        [arf, ix] = sort(arf);
        if arf(1) < fbest, fbest = arf(1);  xbest = arx(:,ix(1)); end
        fhist(end+1) = fbest; %#ok<AGROW>

        if gen == 1 && std(arf) < 1e-12
            warning('dacq_cmaes:flat', ...
                'objective flat in generation 1 -- raise sigma0 (currently %.2f)', sigma0);
        end

        xold  = xmean;
        xmean = arx(:,ix(1:mu)) * w;

        ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold)/sigma;
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*nfe/lambda))/chiN < 1.4 + 2/(n+1);
        pc = (1-cc)*pc + hsig*sqrt(cc*(2-cc)*mueff) * (xmean-xold)/sigma;

        artmp = (arx(:,ix(1:mu)) - xold*ones(1,mu)) / sigma;
        C = (1-c1-cmu)*C + c1*(pc*pc.' + (1-hsig)*cc*(2-cc)*C) ...
          + cmu * artmp * diag(w) * artmp.';

        sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));

        C = triu(C) + triu(C,1).';
        [Bv, Dv] = eig(C);
        dd = real(diag(Dv));  dd(dd < 1e-20) = 1e-20;
        B = real(Bv);  D = sqrt(dd);
        invsqrtC = B * diag(1./D) * B.';

        if sigma < 1e-9 || max(D)/min(D) > 1e10, break; end
    end
    out.fhist = fhist;  out.nfe = nfe;  out.gen = gen;
end
