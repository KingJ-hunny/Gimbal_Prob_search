%==========================================================================
% SGP4
%   SGP4 progagator
%
%   Given:
%     Tle                     structure     NORAD TLE structure
%     timeFromEpochInMinute     double      Time from epoch (minute)            
%
%   Returned (function value):
%     vecr                double     Propagated position vector (3x1, km)
%     vecrdot             double     Propagated velocity vector (3x1, km/s)
%
%   Related:
%
%       loadTle.m   laod TLE data from textfile.
%      
%   References:
%
%       Felix R. Hoots, Ronald L. Roehrich, Spacetrack report NO.3 (1980)
%
%   This revision:  2011 March 25
%
%==========================================================================

function [vecr, vecrdot] =  sgp4(Tle, timeFromEpochInMinute)
    
    % Propagation constants %
    AE = 1;
    ER = 6378.135;
    ke = 0.0743669161;
    s = 1.01222928018927;
    QOMS2T = 1.88027916e-9;
    k2 = 5.413080e-4;
    k4 = 0.62098875e-6;
    J3 = -0.253881e-5;
    A30 = -J3*AE^3;
    D2R = pi/180;

    % Propagation time %
    tMt0 = timeFromEpochInMinute;

    % TLE load and unit converstion %
    n0 = Tle.meanMotion /1440 * 2*pi;
    e0 = Tle.eccentricity;
    i0 = Tle.inclination * D2R;
    M0 = Tle.meanAnomaly * D2R;
    w0 = Tle.argumentOfPerigee * D2R;
    Omega0 = Tle.ascendingNode * D2R;
    ndot0 = Tle.firstMeanMotionDerivative *2 /1440^2 * 2*pi;
    ndot1 = Tle.secondMeanMotionDerivative *6 /1440^3 * 2*pi;
    B = Tle.dragTerm;
    
    % SGP4 Algorithm %
    a1 =  (ke/n0)^(2/3);                                               
    delta1 = 3/2*k2/a1^2*( 3*(cos(i0))^2-1 ) / (1-e0^2)^(3/2);         
    a0 = a1*(1 - 1/3*delta1 - delta1^2 -134/81*delta1^3 );          
    delta0 = 3/2*k2/a0^2*(3*(cos(i0))^2-1)/(1-e0^2)^(3/2);          
    n0dp = n0/(1+delta0);                                           
    a0dp = a0/(1-delta0);                                           

    perigee = a0dp*(1-e0)-AE;
    isImp = 0;if perigee < ( 220/ER ); isImp = 1; end

    if perigee < 156/ER   
        snew = a0dp*(1-e0) - s + AE;                                
        if perigee < 98/ER       
            snew = 20/ER + AE; 
        end    
        QOMS2T = (QOMS2T^(1/4)+s-snew)^4;
        s = snew;
    end

    theta = cos(i0);                     
    zeta = 1/(a0dp - s);                 
    beta0 = (1-e0^2)^(1/2);
    eta = a0dp*e0*zeta;
    C2 = QOMS2T * zeta^4 * n0dp * ( 1- eta^2 )^(-7/2)...
        * (a0dp * ( 1 + 3/2*eta^2 + 4*e0*eta + e0*eta^3)...
        + 3/2*k2*zeta/(1-eta^2) * (-1/2 + 3/2*theta^2) * (8+24*eta^2+3*eta^4));
    C1 = B*C2;
    C3 = QOMS2T * zeta^5 * A30 * n0dp * AE * sin(i0) / (k2*e0);
    C4 = 2*n0dp*QOMS2T*zeta^4*a0dp*beta0^2*(1-eta^2)^(-7/2)...
        *((2*eta*(1+e0*eta) +1/2*e0 +1/2*eta^3) - 2*k2*zeta/(a0dp*(1-eta^2))...
        *(3*(1-3*theta^2)*(1+3/2*eta^2-2*e0*eta-1/2*e0*eta^3)...
        +3/4*(1-theta^2)*(2*eta^2-e0*eta-e0*eta^3)*cos(2*w0)));
    C5 = 2*QOMS2T*zeta^4*a0dp*beta0^2*(1-eta^2)^(-7/2)...
        *(1+11/4*eta*(eta+e0)+e0*eta^3);


    MDF = M0 + (1 + 3*k2*(-1+3*theta^2)/(2*a0dp^2*beta0^3)...
                  + 3*k2^2*(13-78*theta^2+137*theta^4) / ( 16*a0dp^4*beta0^7 ) )...
                  * n0dp*(tMt0);          
    wDF = w0 + ( -3*k2*(1-5*theta^2)/(2*a0dp^2*beta0^4)...
                 + 3*k2^2*(7-114*theta^2 + 395*theta^4)/(16*a0dp^4*beta0^8)...
                 + 5*k4*(3-36*theta^2 + 49*theta^4 ) / (4*a0dp^4*beta0^8 ))...
                 * n0dp*(tMt0);
    deltaw = B*C3*cos(w0)*(tMt0);   
    deltaM = -2/3*QOMS2T * B * zeta^4 * AE/(e0*eta)...
            * ( (1+eta*cos(MDF))^3 - (1+eta*cos(M0))^3 );

    OmegaDF = Omega0 + ( -3*k2*theta / (a0dp^2*beta0^4)...
                         + 3*k2^2*( 4*theta - 19*theta^3 ) / (2*a0dp^4 * beta0^8)...
                         + 5*k4*theta*(3 - 7*theta^2) / (2*a0dp^4*beta0^8) )...
                         *n0dp*(tMt0);

    Omega = OmegaDF - 21/2*(n0dp*k2*theta)/(a0dp^2*beta0^2)*C1*(tMt0)^2;        

    if ~isImp 
        % Higer than 220 km
        D2 = 4*a0dp*zeta*C1^2;
        D3 = 4/3*a0dp*zeta^2*(17*a0dp + s)*C1^3;
        D4 = 2/3*a0dp*zeta^3*(221*a0dp + 31*s)*C1^4;
        
        Mp = MDF + deltaw + deltaM;
        w = wDF - deltaw - deltaM;
        IL = Mp + w + Omega + n0dp* ( 3/2 * C1 * (tMt0)^2 ...
            + (D2 + 2*C1^2)*(tMt0)^3 ...
            +1/4*(3*D3 + 12*C1*D2 + 10*C1^3)*(tMt0)^4 ...
            +1/5*(3*D4 + 12*C1*D3 + 6*D2^2+30*C1^2*D2+15*C1^4)*(tMt0)^5);       
        e = e0 - B*C4*(tMt0) - B*C5*(sin(Mp) - sin(M0));
        a = a0dp*(1-C1*(tMt0) - D2*(tMt0)^2 - D3*(tMt0)^3 - D4*(tMt0)^4)^2;        
    else
        % Lower than 220 km
        Mp = MDF;
        w = wDF;
        IL = Mp + w + Omega + n0dp* ( 3/2 * C1 * (tMt0)^2);
        e = e0 - B*C4*(tMt0);
        a = a0dp*(1-C1*(tMt0))^2;                
    end   
    
    beta = sqrt(1-e^2);
    n = ke/a^(3/2);

    axN = e*cos(w);
    ILL = A30*sin(i0) / (8*k2*a*beta^2) * (e*cos(w)) * ((3+5*theta)/(1+theta));
    ayNL = A30*sin(i0) / (4*k2*a*beta^2);
    ILT = IL + ILL;
    ayN = e*sin(w) + ayNL;

    % Solving kepler equation
    U = ILT - Omega;
    EPw = U;
    DeltaEPw = (U - ayN*cos(EPw) + axN*sin(EPw) - EPw)...
               / (-ayN*sin(EPw) - axN*cos(EPw) + 1);
    while abs(DeltaEPw) > 1e-6    
        EPw = EPw + DeltaEPw;
        DeltaEPw = (U - ayN*cos(EPw) + axN*sin(EPw) - EPw)...
                  / (-ayN*sin(EPw) - axN*cos(EPw) + 1);             
    end

    ecosE = axN*cos(EPw) + ayN*sin(EPw);
    esinE = axN*sin(EPw) - ayN*cos(EPw);
    eL = ( axN^2 + ayN^2 ) ^ (1/2);
    PL = a*(1-eL^2);
    r = a*(1-ecosE);
    rdot = ke*sqrt(a)/r*esinE;
    rfdot = ke*sqrt(PL)/r;
    cosu = a/r * (cos(EPw) - axN + ayN*esinE / (1 + sqrt(1 - eL^2) ) );
    sinu = a/r * (sin(EPw) - ayN - axN*esinE / (1 + sqrt(1 - eL^2) ) );
    u = atan2(sinu, cosu);

    Deltar = k2/(2*PL)*(1-theta^2)*cos(2*u);
    Deltau = -k2/(4*PL^2)*(7*theta^2 - 1) * sin(2*u);
    DeltaOmega = 3*k2*theta/(2*PL^2) * sin(2*u);
    Deltai = 3*k2*theta/(2*PL^2)*sin(i0)*cos(2*u);
    Deltardot = -k2*n/PL*(1-theta^2)*sin(2*u);
    Deltarfdot = k2*n/PL* ( (1-theta^2)*cos(2*u) - 3/2*(1-3*theta^2) );

    rk = r*(1-3/2*k2*sqrt(1-eL^2)/PL^2*(3*theta^2-1)) + Deltar;
    uk = u + Deltau;
    Omegak = Omega + DeltaOmega;
    ik = i0 + Deltai;
    rdotk = rdot + Deltardot;
    rfdotk = rfdot + Deltarfdot;

    vecM = [-sin(Omegak)*cos(ik);
            cos(Omegak)*cos(ik);
            sin(ik)];
    vecN = [cos(Omegak);
            sin(Omegak);
            0];    

    vecU = vecM*sin(uk) + vecN*cos(uk);
    vecV = vecM*cos(uk) - vecN*sin(uk);

    vecr = rk*vecU * ER / AE;
    vecrdot = (rdotk*vecU + rfdotk * vecV) * ER / AE / 60;

end