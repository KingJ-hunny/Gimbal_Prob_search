function [az, el, rng] = eciToAzElVec(rEci, gmst, station)
%==========================================================================
% eciToAzElVec  Vectorised topocentric azimuth/elevation for MANY ECI points
%               at a single instant (same GMST). Same geometry as
%               computeAzAlt.m (Yoon pipeline: ECI->ECEF->ENU), but for a
%               3xN batch of points.
%
%   Given:
%     rEci     3xN  ECI positions [km]
%     gmst     scalar GMST [rad]
%     station  .latDeg .lonDeg .rEcef (3x1 [km])
%
%   Returned (each 1xN):
%     az   azimuth  [rad] in [0,2pi)
%     el   elevation [rad]
%     rng  slant range [km]
%==========================================================================
    R      = eci2ecefMatrix(gmst);          % ECI->ECEF
    ecef   = R * rEci;                       % 3xN
    los    = ecef - station.rEcef(:);        % 3xN (broadcast subtract)

    lat = deg2rad(station.latDeg);
    lon = deg2rad(station.lonDeg);
    Recef2enu = [ -sin(lon),            cos(lon),           0;
                  -sin(lat)*cos(lon),  -sin(lat)*sin(lon),  cos(lat);
                   cos(lat)*cos(lon),   cos(lat)*sin(lon),  sin(lat) ];
    enu = Recef2enu * los;                   % 3xN

    E = enu(1,:); N = enu(2,:); U = enu(3,:);
    rng = sqrt(sum(los.^2, 1));
    az  = mod(atan2(E, N), 2*pi);
    el  = asin(U ./ rng);
end
