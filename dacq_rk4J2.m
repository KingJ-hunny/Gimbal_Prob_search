function [r, v] = dacq_rk4J2(r, v, h, nsteps)
%==========================================================================
% dacq_rk4J2  Fixed-step RK4 propagation of many states under two-body + J2.
%
%   Vectorised over columns: r,v are 3xM and every column is advanced by the
%   same nsteps of size h. Using one integrator for BOTH the sample cloud and
%   its median means the integration error is common-mode and cancels in the
%   differential (cloud-shape) sense, which is what the FOU depends on.
%
%   Given:  r,v  3xM state [km],[km/s];  h step [s];  nsteps count
%   Return: r,v  3xM state after nsteps*h seconds
%==========================================================================
    for s = 1:nsteps
        [k1r, k1v] = deriv(r,            v           );
        [k2r, k2v] = deriv(r + 0.5*h*k1r, v + 0.5*h*k1v);
        [k3r, k3v] = deriv(r + 0.5*h*k2r, v + 0.5*h*k2v);
        [k4r, k4v] = deriv(r +     h*k3r, v +     h*k3v);
        r = r + (h/6)*(k1r + 2*k2r + 2*k3r + k4r);
        v = v + (h/6)*(k1v + 2*k2v + 2*k3v + k4v);
    end
end

% =========================================================================
function [dr, dv] = deriv(r, v)
% Two-body + J2 acceleration, vectorised over columns.
    mu = 398600.4418; J2 = 1.08262668e-3; Re = 6378.137;
    x = r(1,:); y = r(2,:); z = r(3,:);
    r2 = x.*x + y.*y + z.*z;
    rm = sqrt(r2);
    r3 = r2 .* rm;
    r5 = r3 .* r2;

    c  = -1.5 * J2 * mu * Re^2 ./ r5;
    zr = 5 * (z.*z) ./ r2;

    dv = [ -mu*x./r3 + c.*x.*(1 - zr);
           -mu*y./r3 + c.*y.*(1 - zr);
           -mu*z./r3 + c.*z.*(3 - zr) ];
    dr = v;
end
