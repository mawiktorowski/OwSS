close all;
clear all;
format long e;

%param;

[param, ~] = parametry();

thetaRad = 230 * 2 * pi / 360;
dVE = 3;
u = [0; 0; 0; 0; 0; 0];
%u = [0.2; 0.2; 0.2; 0; 0; 0];

zd = [thetaRad; dVE; u];

T = 1;
tau = [0 T/3 2*T/3 T];
ep=1e-6;
rho = 1;

var = obliczenia(param.h0, tau);

[ dQ, grad ] = mexSprawdzenieGrad(zd, param, var, ep, rho);

disp('Porównanie');
disp ([dQ, grad]);
disp('Procentowo (%)');
disp (abs((dQ - grad) ./ grad) * 100);
