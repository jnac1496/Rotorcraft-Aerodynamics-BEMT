function AxialFResults = AxialFlight(y, c, R, thetaC, thetaT, dy, Nb, Omega, vc, Treq, rho, Tmp, swirlFlag)
% INPUTS
% y : the span coordinates from the root cut-out for all the elements
% midpoints (vector of size 1xn) [m]
% c : the chord at distribution at every point of y (vector of size 1xn) [m]
% R : disc radius (scalar) [m]
% thetaC : Collective angle (scalar) [rad]
% thetaT : twist angle distribution at each y location (vector size 1xn) [rad]
% dy : span length of element of the blade (vector of size 1xn) [m]
% Nb : number of blades (scalar)
% Omega: rotor angular speed (scalar) [rad/s]
% vc : climb velocity (scalar) [m]
% Treq: thrust required from the rotor (scalar) [N]
% rho : air density (scalar) [kg/m^3]
% Tmp : air temperature (scalar) [C]
% swirlFlag : boolean variable to include or not swirl effect

% Swirl is by default deactivated
if nargin < 13
    swirlFlag = false;
end

% Global quantities
Vtip = Omega*R;
a = sqrt(1.4*287*(Tmp+273.15));

% Non-dimensional quantities
r = y/R;
sigmaLocal = Nb*c/(pi*R);

% Initial inflow velocity estimate
cla = 2*pi;
theta = thetaC + thetaT;
lambdac = vc/Vtip;
lambda = -(sigmaLocal.*cla/16 - lambdac/2) + ...
    (sqrt((sigmaLocal.*cla/16 - lambdac/2).^2 + sigmaLocal.*cla/8.*theta.*r));
v = lambda*Vtip; % Total flow velocity through disc (in z direction)

% Iteration parameters
n = length(y);
epsilon = 1e-6;
iter_max = 1000;
rf = 0.03;
vmin = 1e-6;
Tmin = 1e-12;

% Preallocation
cl = zeros(1,n);
cd = zeros(1,n);
phi = zeros(1,n);
alpha = zeros(1,n);
F = zeros(1,n);
f = zeros(1,n);
ap = zeros(1,n);
u = zeros(1,n); % tangential velocity or velocity in x direction
Ma = zeros(1,n);
dT = zeros(1,n);
rel_change = zeros(iter_max,n);
res = ones(iter_max,n);
iter = ones(1,n);

% BEMT iteration
for i = 1:n
    while iter(i) <= iter_max
        % Swirl correction
        if ~swirlFlag
            ap(i) = 0;
        else
            swirlRatio = v(i)*(v(i)-vc)/(Omega*y(i))^2;
            determinant = 1 - 4*swirlRatio;
            if determinant < 0
                warning(['Swirl model invalid at element %d. ', ...
                         'Setting swirl correction to zero.'], i);
                ap(i)=0;
            else
                ap(i)=0.5*(1-sqrt(determinant));
            end
        end
        % Relative velocity
        u(i) = Omega*y(i)*(1-ap(i));
        U2 = u(i)^2 + v(i)^2; %Velocity magnitude squared
        % Inflow angle
        phi(i) = atan2(v(i),u(i));
        % Prandtl tip loss
        if r(i) >= 1
            error('Blade element located at r >= 1.');
        end
        f(i) = Nb/2*((1-r(i))/(sin(phi(i))*r(i)));
        F(i) = (2/pi)*acos(exp(-f(i)));
        % Airfoil coefficients
        alpha(i) = theta(i)-phi(i);
        Ma(i) = sqrt(U2)/a;
        [cl(i),cd(i),~] = cpcrcm(alpha(i),Ma(i));
        % Blade element thrust
        dTdybe = 0.5*Nb*rho*c(i)*U2*...
            (cos(phi(i))*cl(i)-sin(phi(i))*cd(i));
        % Momentum thrust
        dTdym = 4*rho*pi*v(i).*(v(i)-vc).*y(i)*F(i);
        % Velocity update
        v_new = vc/2 + sqrt(dTdybe/(4*pi*rho*y(i)*F(i)) + vc^2/4);
        v_next = v(i)+rf*(v_new-v(i));
        % Convergence
        rel_change(iter(i),i) = ...
            abs(v_next-v(i))/max(abs(v_next),vmin);
        res(iter(i),i) = ...
            abs(dTdybe-dTdym)/max(abs(dTdym),Tmin);
        if rel_change(iter(i),i) < epsilon && ...
           res(iter(i),i) < 10*epsilon
            v(i) = v_next;
            break
        end
        v(i)=v_next;
        iter(i)=iter(i)+1;
        if iter(i)>iter_max
            fprintf(['WARNING: no convergence at section %d ', ...
                '(res=%.3e, rel=%.3e)\n'],...
                i,...
                res(iter_max,i),...
                rel_change(iter_max,i));
        end
    end
end

% Recalculate final aerodynamic state using converged vi
lambda = v/Vtip;
lambdai = lambda - lambdac;
vi = lambdai*Vtip;
for i = 1:n
    if ~swirlFlag
        ap(i) = 0;
    else    
        swirlRatio = v(i)*(v(i)-vc)/(Omega*y(i))^2;
        determinant = 1 - 4*swirlRatio;
        if determinant < 0
            warning(['Swirl model invalid at element %d. ', ...
                     'Setting swirl correction to zero.'], i);
            ap(i)=0;
        else
            ap(i)=0.5*(1-sqrt(determinant));
        end 
    end
    u(i) = Omega*y(i)*(1-ap(i));
    U2 = u(i)^2 + v(i)^2;
    phi(i)=atan2(v(i),u(i));
    f(i)=Nb/2*((1-r(i))/(sin(phi(i))*r(i)));
    F(i)=(2/pi)*acos(exp(-f(i)));
    alpha(i)=theta(i)-phi(i);
    Ma(i)=sqrt(U2)/a;
    [cl(i),cd(i),~]=cpcrcm(alpha(i),Ma(i));
    dTdybe = 0.5*Nb*rho*c(i)*U2*...
        (cos(phi(i))*cl(i)-sin(phi(i))*cd(i));
    dT(i)=dTdybe*dy(i);
end

% Thrust calculation
T=sum(dT);
dcT=dT/(rho*Omega^2*pi*R^4);
cT=sum(dcT);
% Thrust difference
if Treq > 0
    DeltaT=T-Treq;
else
    DeltaT=0;
end

% Torque calculation
dQ = 0.5*Nb*rho*c.*(u.^2+v.^2).*(sin(phi).*cl + cos(phi).*cd).*y.*dy;
Q = sum(dQ);
cQ = Q/(rho*Vtip^2*pi*R^3);

% Power calculation
% Just for checking
dPL_old = 0.5*Nb*rho*Omega*c.*(u.^2+v.^2).*sin(phi).*cl.*y.*dy;
PL_old = sum(dPL_old);
cPL_old = PL_old/(rho*Vtip^3*pi*R^2);
fprintf('Lift power coeff (direct): %.6e\n', cPL_old);

% New Generealized deinfition
dTL = Nb*(0.5*rho*c.*(u.^2+v.^2).*cl.*dy).*cos(phi); % Lift component for thrust
% Induced power contribution
dPi = dTL.*vi; 
Pi = sum(dPi);
cPi = Pi/(rho*Vtip^3*pi*R^2);
% Climb power contribution
dPc = dTL.*vc; 
Pc = sum(dPc);
cPc = Pc/(rho*Vtip^3*pi*R^2);

% Swirl effect contribution
dPs = ap.*dTL.*v/(1-ap);
Ps = sum(dPs);
cPs = Ps/(rho*Vtip^3*pi*R^2);

% Total Lift power contribution to shaft power required
dPL = dPi + dPc + dPs;
PL = sum(dPL);
cPL = PL/(rho*Vtip^3*pi*R^2);

% Profile power contribution to shaft power required
dPo = 0.5*Nb*rho*Omega*c.*(u.^2+v.^2).*cos(phi).*cd.*y.*dy;
Po = sum(dPo);
cPo = Po/(rho*Vtip^3*pi*R^2);

% Total power
dP = dPL + dPo;
P = sum(dP);
cP = P/(rho*Vtip^3*pi*R^2);

lambdai_ideal = -lambdac/2 + sqrt(lambdac^2/4 + cT/2);
cPIdeal  = cT*(lambdac + lambdai_ideal); % Generalized ideal cP
% If lambdac = 0 then it reduces to cT^(3/2)/sqrt(2)
cPiIdeal = cT*lambdai_ideal; % -> cT^(3/2)/sqrt(2) at lambdac = 0
% Defined to be meaningful during climb otherwise k<1
k = cPi/cPiIdeal;

% Performance
if vc == 0
    % Figure of Merit for hover
    FM = cPIdeal/(cP);
    etaC = NaN;
else
    % Climb efficiency
    etaC = cPIdeal/cP;
    FM = NaN;
end


% Output
AxialFResults = struct(...
    'InflowVelocity',v,...
    'InflowCoeff',lambda,...
    'InducedVelocity', vi,...
    'InducedInflowCoeff',lambdai,...
    'Thrust',T,...
    'ThrustCoeff',cT,...
    'ThrustDifference',DeltaT,...
    'ThrustDistri', dT,...
    'Torque', Q,...
    'TorqueCoeff', cQ,...
    'TorqueDistri', dQ,...
    'InducedPower', Pi,...
    'InducedPowerCoeff', cPi,...
    'ClimbPower', Pc,...
    'ClimbPowerCoeff', cPc,...
    'SwirlPower', Ps,...
    'SwirlPowerCoeff', cPs,...
    'LiftPower', PL,...
    'LiftPowerCoeff', cPL,...    
    'ProfilePower', Po,...
    'ProfilePowerCoeff', cPo,...
    'TotalPower', P,...
    'TotalPowerCoeff', cP,...
    'FigureOfMerit', FM,...
    'ClimbEfficiency', etaC,...
    'InducedPowerK', k,...
    'SwirlEnabled', swirlFlag,...
    'SwirlCoeff',ap);
end