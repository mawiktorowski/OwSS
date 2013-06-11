function out = x0(LEO, theta, alpha, dV, m)
% wyliczanie warunkow poczatkowych

global mu restmu aE C2 D

r = (aE + LEO) / D;
thetaRad = theta * 2 * pi / 360;
alphaRad = alpha * 2 * pi / 360;
V = sqrt(restmu/r);

out = zeros(1,5);
out(1) = r * cos(thetaRad) - mu;
out(2) = r * sin(thetaRad);
out(3) = dV * sin(alphaRad - thetaRad) - V * sin(thetaRad) + out(2);
out(4) = dV * cos(alphaRad - thetaRad) + V * cos(thetaRad) - out(1);
out(5) = m / 1000 * exp(- C2 * dV);
end