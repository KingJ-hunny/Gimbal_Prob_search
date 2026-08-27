function [R, Rdot] = eci2ecefMatrix(gmst)
%==========================================================================
% eci2ecefMatrix  Rotation from ECI (TEME) to ECEF about the pole by GMST.
%                 r_ecef = R * r_eci.  This is the simple Earth-rotation
%                 model used in the Yoon tracking derivation (ECEF_T_ECI),
%                 neglecting precession/nutation/polar motion.
%
%   Given:  gmst   Greenwich Mean Sidereal Time [rad]
%   Return: R      3x3 rotation matrix (ECI -> ECEF)
%           Rdot   3x3 time derivative of R [1/s] (uses Earth rotation rate)
%==========================================================================
    c = cos(gmst); s = sin(gmst);
    R = [  c,  s, 0;
          -s,  c, 0;
           0,  0, 1 ];

    if nargout > 1
        we = 7.292115146706979e-5;   % Earth rotation rate [rad/s]
        Rdot = we * [ -s,  c, 0;
                      -c, -s, 0;
                       0,  0, 0 ];
    end
end
