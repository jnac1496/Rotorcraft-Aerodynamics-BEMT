function ThrustIter = AxialThrustIter(y, c, R, thetaT, dy, Nb, Omega, vc, Treq, rho, Tmp, swirlFlag)
% INPUTS
% y : the span coordinates from the root cut-out for all the elements ends and its's 
% midpoints (vector of size 1x2*n+1) [m]
% c : the chord at distribution at every point of y (vector of size
% 1x2*n+1) [m]
% R : disc radius (scalar) [m]
% thetaT : twist angle distribution at each y location (vector size 1x2*n+1) [deg]
% dy : span length of element of the blade (vector of size 1xn) [m]
% Nb : number of blades (scalar)
% vc : climb velocity (scalar) [m]
% Omega: rotor angular speed (scalar) [rad/s]
% Treq: thrust required from the rotor (scalar) [N]
% rho : air density (scalar) [kg/m^3]
% Tmp : air temperature (scalar) [C]


% OUTPUTS

if Treq <= 0
    error('Treq must be greater than zero for thrust trim.');
end
if nargin < 12
    swirlFlag = false;
end

% Angles in deg to rad
% provided
thetaT = deg2rad(thetaT);

% Relevant global variables
Vtip = Omega*R;
cla = 2*pi;
%a = sqrt(1.4*287*(Tmp+273.15));

% Non-dimensional quantities
r = y/R; % vector (1x2*n+1)
sigmaLocal = Nb*c/(pi*R); % vector (1x2*n+1)

% Collective angle estimation routine
% Uniform induced velocity guess from Momentum Theory (0th)
vi_0 = sqrt(Treq/(2*rho*pi*R^2));
% Collective angle estimation from uniform vi(0th)
kn1 = trapz(y,thetaT.*y.^2.*c);
kn2 = trapz(y,vi_0/Omega*y.*c);
kden = trapz(y, y.^2.*c);
thetaC = Treq/(0.5*Nb*rho*Omega^2*cla*kden) - kn1/kden + kn2/kden;
theta = thetaC + thetaT; % vector (1x2*n+1)
% Induced velocity estimation from Angular Momentum Th (1st)
lambdai = (sigmaLocal.*cla/16).*(-1 + sqrt(1 + ((32.*theta.*r)./(sigmaLocal.*cla)))); % vector (1x2*n+1)
vi = lambdai*Vtip; % vector (1x2*n+1)
% Collective angle estimation from variable vi (1st)
kn2 = trapz(y,vi./Omega.*y.*c);
thetaC = Treq/(0.5*Nb*rho*Omega^2*cla*kden) - kn1/kden + kn2/kden;


% Prepare data for Hover thrust calculation function
ybem = y(2:2:end); %vector (1xn)
cbem = c(2:2:end); %vector (1xn)
thetaTbem = thetaT(2:2:end); %vector (1xn)

% Tolerance
epsilon = 1e-4;
% First hover thrust calculation
AxialFResults = AxialFlight(ybem, cbem, R, thetaC, thetaTbem, dy, Nb, Omega, vc, Treq, rho, Tmp, swirlFlag);
% Collective angle search for required thrust
if abs(AxialFResults.ThrustDifference)/Treq <= epsilon
    fprintf('Thrust required reached with collective angle = %.4f deg\n', rad2deg(thetaC));
else
    % Bracket search
    dtheta = deg2rad(1);  % initial search step
    iterBracketMax = 50;
    iterBracket = 1;
    bracketFound = false;
    if AxialFResults.ThrustDifference < 0
        thetaC_low = thetaC;
        thetaC_test = thetaC + dtheta;
        while iterBracket <= iterBracketMax
            AxialFResults = AxialFlight(ybem, cbem, R, thetaC_test, thetaTbem, dy, Nb, Omega, vc, Treq, rho, Tmp, swirlFlag);
            if AxialFResults.ThrustDifference >= 0
                thetaC_high = thetaC_test;
                bracketFound = true;
                break
            end
            thetaC_low = thetaC_test;
            dtheta = 2*dtheta;
            thetaC_test = thetaC_test + dtheta;
            iterBracket = iterBracket + 1; % Increment iteration counter
        end
    else 
        thetaC_high = thetaC;
        thetaC_test = thetaC - dtheta;
        while iterBracket <= iterBracketMax
            AxialFResults = AxialFlight(ybem, cbem, R, thetaC_test, thetaTbem, dy, Nb, Omega, vc, Treq, rho, Tmp, swirlFlag);
            if AxialFResults.ThrustDifference <= 0
                thetaC_low = thetaC_test;
                bracketFound = true;
                break
            end
            thetaC_high = thetaC_test;
            dtheta = 2*dtheta;
            thetaC_test = thetaC_test - dtheta;
            iterBracket = iterBracket + 1; % Increment iteration counter          
        end    
    end
    % Bisection method
    iterBisectMax = 100;
    iterBisect = 1;
    if ~bracketFound
        error('Unable to bracket required thrust.');
    else
        while iterBisect <= iterBisectMax
            thetaC = 0.5*(thetaC_low + thetaC_high);
            AxialFResults = AxialFlight(ybem, cbem, R, thetaC, thetaTbem, dy, Nb, Omega, vc, Treq, rho, Tmp, swirlFlag);
            thrust_error = abs(AxialFResults.ThrustDifference)/Treq;
            if thrust_error <= epsilon        
                fprintf('Thrust required reached with collective angle = %.4f deg\n', rad2deg(thetaC));
                break;
            elseif AxialFResults.ThrustDifference < 0
                thetaC_low = thetaC;
            else
                thetaC_high = thetaC;
            end
            iterBisect = iterBisect + 1; % Increment iteration counter
            if iterBisect > iterBisectMax
                error('Bisection method failed to converge after %d iterations.', ...
                    iterBisectMax);
            end
        end
    end
end
ThrustIter.thetaC = thetaC;
ThrustIter.ThrustDifference = AxialFResults.ThrustDifference;
ThrustIter.hoverResults = AxialFResults;
ThrustIter.RelativeThrustError = abs(AxialFResults.ThrustDifference)/Treq;
end