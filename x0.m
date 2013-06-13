function out = x0(LEO, theta, alpha, dV, m)
% wyliczanie warunkow poczatkowych

% parametry które nie s¹ konieczne
%LEO m

% parametry konieczne do istnienia modelu (zmiennr decyzyjne)
%theta alpha dV

global mu C2 rE VE

thetaRad = theta * 2 * pi / 360;
alphaRad = alpha * 2 * pi / 360;

out = zeros(1,5);
out(1) = rE * cos(thetaRad) - mu;
out(2) = rE * sin(thetaRad);
out(3) = dV * sin(alphaRad - thetaRad) - VE * sin(thetaRad) + out(2);
out(4) = dV * cos(alphaRad - thetaRad) + VE * cos(thetaRad) - out(1);
out(5) = m * exp(- C2 * dV);
end