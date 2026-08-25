function geometry = BladeGeom(ro, R, co, lambda, Dtheta, Nb, q, p, n)
% INPUTS
% ro : blade cut-put (scalar) [m]
% R : disc radius (scalar) [m]
% co : root chord (scalar) [m]
% lambda : tapper ratio (scalar)
% Dtheta : Delta theta between tip and root pitch (negative scalar) [deg]
% Nb : number of blades (scalar)
% q : chord distribution 5th order Bezier curve control points (vector ascendin
% order of size 1x4)
% p : twist distribution 5th order Bezier curve control points (vector ascendin
% order of size 1x4)
% n : number of blade elements for BEM

% OUTPUTS
% geometry.y : the span coordinates from the root cut-out for all the elements 
% and its's midpoints (vector of size 1x2*n+1) [m]
% geometry.c : the chord at distribution at every point of y (vector of size
% 1x2*n+1) [m]
% geometry.thetaTwist : twist distribution at every point of y (vector of size
% 1x2*n+1) [deg]
% geometry.yBem : span coordinates of mid-points of each blade element (vector 
% of size 1xn) [m]
% geometry.cBem : chords at every yBem location (vector of size 1xn) [m]
% geometry.thetaTwistBem : twist angle at every yBem location (vector of 
% size 1xn) [deg]
% geometry.dy : span lenght of each of the n elements (vector of size 1xn) [m]
% geometry.sigmaLocalBem : local solidity at every yBem location (vector of 
% size 1xn)
% geometry.sigmaArea : solidity of the aerodynamic surfaces of the disc
% (scalar)
% geometry.cMean : mean chord of the blade [m]
% geometry.sigmaNominal : nominal solidity considering cMean

% Chord distribution
%dy = (R-ro)/n; 
%ye =linspace(dy+ro,R,n); %
y = linspace(ro,R,2*n+1); %ro:dy/2:R;

s = (y -ro)./(R - ro); % Bezier function variable

% Implement Bezier curve for distribitions
q = [0 q 1];
p = [0 p 1];
Bord = length(q)-1; % order of Bezier function 
c = zeros(1,length(s));
theta = zeros(1,length(s));

for i = 1:length(s)
    hc = 0;
    hp = 0;
    for j = 0:Bord
        hc = hc + q(j+1)*nchoosek(Bord,j)*s(i)^j*(1-s(i))^(Bord-j);
        hp = hp + p(j+1)*nchoosek(Bord,j)*s(i)^j*(1-s(i))^(Bord-j);
    end
    c(i) = co*(1+(lambda-1)*hc);
    theta(i) = Dtheta*hp;
end

% Mid element coordinates, chord and twist for BEM
ybem = y(2:2:end);
cbem = c(2:2:end);
thetabem = theta(2:2:end);
% Span of each element along the blade
yEdges = y(1:2:end);
dy = diff(yEdges);

% Compute solidity
% Solidity for aerodynamic surfaces
sigmaArea = Nb*trapz(y,c)/(pi*R^2);

% Mean chord over the active blade span
cMean = trapz(y,c)/(R-ro);

% Nominal rotor solidity
sigmaNominal = Nb*cMean/(pi*R);

% Local BEM solidity
sigmaLocalBem = Nb*cbem/(pi*R);

% Output
geometry.y = y;
geometry.c = c;
geometry.thetaTwist = theta;
geometry.yBem = ybem;
geometry.cBem = cbem;
geometry.thetaTwistBem = thetabem;
geometry.dy = dy;
geometry.sigmaLocalBem = sigmaLocalBem;
geometry.sigmaArea = sigmaArea;
geometry.cMean = cMean;
geometry.sigmaNominal = sigmaNominal;
end