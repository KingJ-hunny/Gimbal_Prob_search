function fou = projectFoU(SigmaEci, rEci, jd, station, opt)
%==========================================================================
% projectFoU  Project a 3-D ECI position covariance onto the plane of sky
%             (the local az/alt tangent plane) to obtain the angular
%             Field-of-Uncertainty (FOU) seen from a ground station.
%
%   Pipeline (Kim et al., Optics Communications 620, 2026, Sec. 4.1):
%     1. Rotate the ECI covariance into the station ENU frame.
%     2. Keep only the two components perpendicular to the line of sight
%        (radial/range error does not move the pointing direction):
%           x-axis = cross-elevation  (increasing azimuth, on-sky)
%           y-axis = elevation        (increasing altitude)
%     3. Divide the linear 2-D covariance by range^2 -> angular covariance.
%        The orbit term therefore subtends ~ sigma_pos / range, so the FOU
%        SHRINKS with range and GROWS as the pass climbs toward zenith.
%     4. Optionally add an isotropic angular term (attitude / boresight /
%        thermoelastic) in quadrature: Sigma_ang += sigma_iso^2 * I.
%
%   Given:
%     SigmaEci  3x3xN  ECI position covariance [km^2] (from generateCovariance)
%     rEci      3xN    satellite ECI position [km]
%     jd        1xN    absolute Julian Date of each sample (for GMST)
%     station   struct  .latDeg .lonDeg .altKm   (ground station)
%     opt       struct (optional):
%                 .nSigma       ellipse sigma-level (default 3)
%                 .sigmaIsoRad  isotropic angular 1-sigma [rad] (default 0)
%                 .nEllipse     boundary sample count (default 120)
%
%   Returned struct fou (all 1xN unless noted):
%     .azDeg .altDeg .rangeKm     pointing + slant range
%     .SigmaAng   2x2xN  angular covariance [rad^2] in [crossEl, El] basis
%     .aRad .bRad         nSigma semi-major / semi-minor angular half-axes [rad]
%     .k                  axis ratio b/a  (elliptical FOU anisotropy)
%     .angleDeg           major-axis orientation, from the cross-elevation
%                         (azimuth) axis toward the elevation axis [deg]
%     .ellipAzDeg .ellipAltDeg   (nEllipse x N) ellipse boundary drawn in the
%                         az/alt plane (deg), ready for plotting on a sky map
%==========================================================================

    if nargin < 5 || isempty(opt), opt = struct(); end
    if ~isfield(opt,'nSigma')      || isempty(opt.nSigma),      opt.nSigma = 3;   end
    if ~isfield(opt,'sigmaIsoRad') || isempty(opt.sigmaIsoRad), opt.sigmaIsoRad = 0; end
    if ~isfield(opt,'nEllipse')    || isempty(opt.nEllipse),    opt.nEllipse = 120; end

    if ~isfield(station,'rEcef') || isempty(station.rEcef)
        station.rEcef = geodetic2ecef(station.latDeg, station.lonDeg, station.altKm);
    end

    N = size(rEci, 2);
    fou.azDeg    = zeros(1, N);
    fou.altDeg   = zeros(1, N);
    fou.rangeKm  = zeros(1, N);
    fou.SigmaAng = zeros(2, 2, N);
    fou.aRad     = zeros(1, N);
    fou.bRad     = zeros(1, N);
    fou.k        = zeros(1, N);
    fou.angleDeg = zeros(1, N);
    fou.ellipAzDeg  = zeros(opt.nEllipse, N);
    fou.ellipAltDeg = zeros(opt.nEllipse, N);

    th = linspace(0, 2*pi, opt.nEllipse);

    for k = 1:N
        gmst = gmstRad(jd(k));
        [az, alt, rng, ~, Reci2enu] = computeAzAlt(rEci(:,k), gmst, station);

        % ECI covariance -> ENU
        Senu = Reci2enu * SigmaEci(:,:,k) * Reci2enu.';

        % Sky-tangent basis (perpendicular to LOS), expressed in ENU:
        %   eX : cross-elevation (increasing azimuth)
        %   eY : elevation       (increasing altitude)
        eX = [ cos(az);            -sin(az);            0        ];
        eY = [-sin(alt)*sin(az);   -sin(alt)*cos(az);   cos(alt) ];
        P  = [eX.'; eY.'];                       % 2x3

        Slin = P * Senu * P.';                   % linear 2-D covariance [km^2]
        Sang = Slin / rng^2;                     % angular covariance [rad^2]
        Sang = 0.5*(Sang + Sang.');              % symmetrize
        Sang = Sang + opt.sigmaIsoRad^2 * eye(2);% add isotropic angular term

        % Eigen-decomposition -> ellipse axes/orientation
        [V, D] = eig(Sang);
        lam = diag(D);
        [lam, idx] = sort(lam, 'descend');
        V = V(:, idx);
        lam = max(lam, 0);

        aRad = opt.nSigma * sqrt(lam(1));        % semi-major (angular)
        bRad = opt.nSigma * sqrt(lam(2));        % semi-minor (angular)
        vmaj = V(:,1);

        % Ellipse boundary in [crossEl, El] angular coords
        bx = aRad*cos(th)*vmaj(1) + bRad*sin(th)*(-vmaj(2));
        by = aRad*cos(th)*vmaj(2) + bRad*sin(th)*( vmaj(1));

        % Map angular offsets to az/alt degrees for a sky-plane plot:
        %   d(alt) = by ;  d(az) = bx / cos(alt)   (on-sky cross-El = cosalt*dAz)
        fou.ellipAzDeg(:,k)  = rad2deg(az)  + rad2deg(bx / cos(alt));
        fou.ellipAltDeg(:,k) = rad2deg(alt) + rad2deg(by);

        fou.azDeg(k)    = rad2deg(az);
        fou.altDeg(k)   = rad2deg(alt);
        fou.rangeKm(k)  = rng;
        fou.SigmaAng(:,:,k) = Sang;
        fou.aRad(k)     = aRad;
        fou.bRad(k)     = bRad;
        fou.k(k)        = bRad / aRad;
        fou.angleDeg(k) = rad2deg(atan2(vmaj(2), vmaj(1)));
    end
end
