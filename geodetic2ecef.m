function rEcef = geodetic2ecef(latDeg, lonDeg, altKm)
%==========================================================================
% geodetic2ecef  WGS-84 geodetic (lat, lon, alt) -> ECEF position [km].
%
%   Given:  latDeg  geodetic latitude  [deg]
%           lonDeg  longitude          [deg, east positive]
%           altKm   height above WGS-84 ellipsoid [km]
%   Return: rEcef   ECEF position vector (3x1) [km]
%==========================================================================
    a  = 6378.137;                 % WGS-84 semi-major axis [km]
    f  = 1/298.257223563;          % flattening
    e2 = f*(2 - f);                % first eccentricity squared

    lat = deg2rad(latDeg);
    lon = deg2rad(lonDeg);

    N = a / sqrt(1 - e2*sin(lat)^2);          % prime-vertical radius of curvature
    rEcef = [ (N + altKm)*cos(lat)*cos(lon);
              (N + altKm)*cos(lat)*sin(lon);
              (N*(1 - e2) + altKm)*sin(lat) ];
end
