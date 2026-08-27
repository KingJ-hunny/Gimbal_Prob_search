function Tle = parseTLE(line1, line2)
%==========================================================================
% parseTLE  Parse a two-line element set (TLE) into the struct expected by
%           sgp4.m, plus absolute epoch information.
%
%   Given:
%     line1   char/string   TLE line 1 (69 columns, standard format)
%     line2   char/string   TLE line 2 (69 columns, standard format)
%
%   Returned:
%     Tle     struct with fields consumed by sgp4.m:
%               .meanMotion                 [rev/day]
%               .eccentricity               [-]
%               .inclination                [deg]
%               .meanAnomaly                [deg]
%               .argumentOfPerigee          [deg]
%               .ascendingNode              [deg]  (RAAN)
%               .firstMeanMotionDerivative  [rev/day^2]  (raw TLE value = ndot/2)
%               .secondMeanMotionDerivative [rev/day^3]  (raw TLE value = nddot/6)
%               .dragTerm                   [1/ER]       (BSTAR)
%             and extra bookkeeping fields:
%               .satnum        satellite catalog number
%               .epochYear     4-digit year
%               .epochDay      day-of-year (fractional)
%               .epochJD       Julian Date (UTC) of the epoch
%
%   Notes:
%     * Fixed-column parsing per the NORAD/Vallado TLE convention.
%     * BSTAR and nddot use the implied-decimal exponential notation
%       (e.g. " 78601-3" -> 0.78601e-3 ; "-19221-5" -> -0.19221e-5).
%==========================================================================

    l1 = char(line1);
    l2 = char(line2);

    % Pad to 69 columns so fixed-column slicing never overruns.
    if numel(l1) < 69, l1(end+1:69) = ' '; end
    if numel(l2) < 69, l2(end+1:69) = ' '; end

    % ---- Line 1 -------------------------------------------------------
    Tle.satnum  = str2double(l1(3:7));

    yy = str2double(l1(19:20));
    if yy < 57                 % TLE two-digit-year rollover convention
        Tle.epochYear = 2000 + yy;
    else
        Tle.epochYear = 1900 + yy;
    end
    Tle.epochDay = str2double(l1(21:32));

    Tle.firstMeanMotionDerivative  = str2double(strtrim(l1(34:43)));  % ndot/2
    Tle.secondMeanMotionDerivative = parseImpliedExp(l1(45:52));      % nddot/6
    Tle.dragTerm                   = parseImpliedExp(l1(54:61));      % BSTAR

    % ---- Line 2 -------------------------------------------------------
    Tle.inclination       = str2double(l2(9:16));
    Tle.ascendingNode     = str2double(l2(18:25));
    Tle.eccentricity      = str2double(['0.' strtrim(l2(27:33))]);
    Tle.argumentOfPerigee = str2double(l2(35:42));
    Tle.meanAnomaly       = str2double(l2(44:51));
    Tle.meanMotion        = str2double(l2(53:63));

    % ---- Epoch as Julian Date (UTC) -----------------------------------
    % JD at 0h Jan 1 of the epoch year, then add (day-of-year - 1).
    Tle.epochJD = jday(Tle.epochYear, 1, 1, 0, 0, 0) + (Tle.epochDay - 1);
end

% -------------------------------------------------------------------------
function val = parseImpliedExp(field)
% Parse TLE implied-decimal exponential fields such as " 78601-3".
    s = strtrim(field);
    if isempty(s) || all(s == ' ')
        val = 0; return;
    end
    mant = s(1:end-2);          % mantissa digits (with optional sign)
    ex   = s(end-1:end);        % signed one-digit exponent, e.g. '-3'
    sgn  = 1;
    if mant(1) == '-'
        sgn = -1; mant = mant(2:end);
    elseif mant(1) == '+'
        mant = mant(2:end);
    end
    val = sgn * str2double(['0.' mant]) * 10^str2double(ex);
end
