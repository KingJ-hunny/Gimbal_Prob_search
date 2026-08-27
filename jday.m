function jd = jday(yr, mon, day, hr, mn, sec)
%==========================================================================
% jday  Julian Date from a Gregorian calendar date/time (UTC).
%       Valid for years 1900-2100 (Vallado, "Fundamentals of Astrodynamics
%       and Applications").
%
%   Given:  yr, mon, day, hr, mn, sec   (day may be fractional)
%   Return: jd   Julian Date [days]
%==========================================================================
    jd = 367*yr ...
       - floor( 7*(yr + floor((mon + 9)/12)) / 4 ) ...
       + floor( 275*mon/9 ) ...
       + day + 1721013.5 ...
       + ((sec/60 + mn)/60 + hr)/24;
end
