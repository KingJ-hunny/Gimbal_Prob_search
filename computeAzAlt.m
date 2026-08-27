function [az, alt, rangeKm, enuHat, Reci2enu] = computeAzAlt(rSatEci, gmst, station)
%==========================================================================
% computeAzAlt  Topocentric azimuth/altitude of a satellite from a fixed
%               ground station, following the Yoon (AE550) pointing pipeline:
%
%                 LOS in ECI  ->  ECEF (R3 by GMST)  ->  local ENU (gimbal
%                 base frame)  ->  spherical angles (az, alt).
%
%   This is the concrete realization of Yoon's frame chain
%   rgim = gim_T_ECEF * ECEF_T_ECI * rECI, with the gimbal base frame taken
%   as the station East-North-Up frame; the resulting (az, alt) are identical
%   to the closed-form eqs. (8)-(9) of the solution for a standard az-alt mount.
%
%   Given:
%     rSatEci  3x1  satellite position in ECI [km]
%     gmst     scalar GMST [rad] at the observation instant
%     station  struct with fields:
%                .latDeg .lonDeg .altKm        (geodetic location)
%                .rEcef  (3x1, optional cache) (ECEF position [km])
%
%   Returned:
%     az        azimuth  [rad], measured from North toward East, in [0,2pi)
%     alt       altitude/elevation [rad], in [-pi/2, pi/2]
%     rangeKm   slant range station->satellite [km]
%     enuHat    3x1 unit line-of-sight in the ENU frame
%     Reci2enu  3x3 rotation ECI -> ENU (useful for covariance projection)
%==========================================================================

    if ~isfield(station, 'rEcef') || isempty(station.rEcef)
        station.rEcef = geodetic2ecef(station.latDeg, station.lonDeg, station.altKm);
    end

    % --- ECI -> ECEF -------------------------------------------------------
    Reci2ecef = eci2ecefMatrix(gmst);
    rSatEcef  = Reci2ecef * rSatEci(:);

    % --- Line of sight in ECEF --------------------------------------------
    losEcef = rSatEcef - station.rEcef(:);
    rangeKm = norm(losEcef);

    % --- ECEF -> ENU at the station ---------------------------------------
    lat = deg2rad(station.latDeg);
    lon = deg2rad(station.lonDeg);
    Recef2enu = [ -sin(lon),            cos(lon),           0;
                  -sin(lat)*cos(lon),  -sin(lat)*sin(lon),  cos(lat);
                   cos(lat)*cos(lon),   cos(lat)*sin(lon),  sin(lat) ];

    enu    = Recef2enu * losEcef;
    enuHat = enu / rangeKm;

    E = enu(1); N = enu(2); U = enu(3);

    % --- Spherical angles --------------------------------------------------
    az  = mod(atan2(E, N), 2*pi);          % from North, clockwise (toward East)
    alt = asin(U / rangeKm);               % elevation above local horizon

    % Full ECI -> ENU rotation, for projecting an ECI covariance onto the sky.
    Reci2enu = Recef2enu * Reci2ecef;
end
