function gmst = gmstRad(jdUT1)
%==========================================================================
% gmstRad  Greenwich Mean Sidereal Time [rad] from Julian Date (UT1~UTC).
%          IAU-1982 polynomial (Vallado Eq. 3-47). Output wrapped to [0,2pi).
%
%   Given:  jdUT1   Julian Date [days]
%   Return: gmst    GMST [rad], in [0, 2*pi)
%==========================================================================
    Tu = (jdUT1 - 2451545.0) / 36525.0;           % Julian centuries from J2000
    gmstSec = 67310.54841 ...
            + (876600*3600 + 8640184.812866)*Tu ...
            + 0.093104*Tu.^2 ...
            - 6.2e-6*Tu.^3;                        % seconds
    gmstDeg = mod(gmstSec, 86400) / 240;           % 86400 s -> 360 deg  =>  /240 s/deg
    gmst = deg2rad(mod(gmstDeg, 360));
end
