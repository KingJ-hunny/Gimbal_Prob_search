function cov = generateCovariance(trk, sigma)
%==========================================================================
% generateCovariance  (FUNCTION 2)  Along-track-elongated orbit covariance.
%
%   Builds the 3x3 position covariance of the orbit-prediction error at every
%   sample of a track, expressed both in the local RIC (Radial / In-track /
%   Cross-track) frame and rotated into ECI.  The error is anisotropic and
%   dominated by the in-track component, matching the SGP4/TLE error structure
%   reported by Kim et al. (Optics Communications 620, 2026) and Easthope
%   (2015): in-track >> cross-track ~ radial.
%
%   RIC frame convention (right-handed):
%     R_hat = r / |r|                        (radial, outward)
%     C_hat = (r x v) / |r x v|              (cross-track, orbit normal)
%     I_hat = C_hat x R_hat                  (in-track, ~ velocity direction)
%   Rotation RIC -> ECI :  Q = [R_hat  I_hat  C_hat]  (columns), so that
%     Sigma_ECI = Q * Sigma_RIC * Q'.
%
%   Given:
%     trk    struct from generateTrack (needs .rEci, .vEci, 3xN)
%     sigma  one of:
%              - struct with scalar 1-sigma std-devs [km]:
%                    .R  radial      (default 0.17)
%                    .I  in-track    (default 0.80)
%                    .C  cross-track (default 0.17)
%                (giving an in-track-dominated diagonal RIC covariance)
%              - struct additionally carrying .SigmaRIC (3x3) full RIC
%                covariance, used verbatim (overrides R/I/C)
%              - a bare 3x3 matrix, taken as a constant RIC covariance
%
%   Returned struct cov:
%     .SigmaRIC   3x3xN  covariance in the RIC frame [km^2]
%     .SigmaEci   3x3xN  covariance rotated into ECI [km^2]
%     .Rhat/.Ihat/.Chat  3xN  RIC basis vectors in ECI
%     .sigmaRIC   [sR sI sC]  the 1-sigma values used [km]
%     .kOrb       min(sC,sR)/sI   position-error axis ratio (Kim et al.)
%==========================================================================

    r = trk.rEci;  v = trk.vEci;
    N = size(r, 2);

    % ---- resolve the RIC covariance template -------------------------
    if isnumeric(sigma) && isequal(size(sigma), [3 3])
        SigmaRIC0 = sigma;
        sR = sqrt(sigma(1,1)); sI = sqrt(sigma(2,2)); sC = sqrt(sigma(3,3));
    elseif isstruct(sigma) && isfield(sigma, 'SigmaRIC') && ~isempty(sigma.SigmaRIC)
        SigmaRIC0 = sigma.SigmaRIC;
        sR = sqrt(SigmaRIC0(1,1)); sI = sqrt(SigmaRIC0(2,2)); sC = sqrt(SigmaRIC0(3,3));
    else
        if ~isstruct(sigma), sigma = struct(); end
        if ~isfield(sigma,'R') || isempty(sigma.R), sigma.R = 0.17; end
        if ~isfield(sigma,'I') || isempty(sigma.I), sigma.I = 0.80; end
        if ~isfield(sigma,'C') || isempty(sigma.C), sigma.C = 0.17; end
        sR = sigma.R; sI = sigma.I; sC = sigma.C;
        SigmaRIC0 = diag([sR^2, sI^2, sC^2]);
    end

    cov.SigmaRIC = repmat(SigmaRIC0, [1 1 N]);
    cov.SigmaEci = zeros(3, 3, N);
    cov.Rhat = zeros(3, N);
    cov.Ihat = zeros(3, N);
    cov.Chat = zeros(3, N);

    for k = 1:N
        rk = r(:,k); vk = v(:,k);
        Rhat = rk / norm(rk);
        h    = cross(rk, vk);
        Chat = h / norm(h);
        Ihat = cross(Chat, Rhat);           % right-handed, ~ along velocity

        Q = [Rhat, Ihat, Chat];             % RIC -> ECI
        cov.SigmaEci(:,:,k) = Q * SigmaRIC0 * Q.';

        cov.Rhat(:,k) = Rhat;
        cov.Ihat(:,k) = Ihat;
        cov.Chat(:,k) = Chat;
    end

    cov.sigmaRIC = [sR, sI, sC];
    cov.kOrb     = min(sC, sR) / sI;
end
