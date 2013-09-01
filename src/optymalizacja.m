close all;
clear all;
format long e;

[param, opts] = parametry();

thetaRad = 230 * 2 * pi / 360;
dVE = 3;
%u = [0; 0; 0; 0; 0; 0];
u = [0.2; 0.2; 0.2; 0; 0; 0];

zd = [thetaRad; dVE; u];

T = 1;
tau = [0 T/3 2*T/3 T];
rho = 1;

var = obliczenia(param.h0, tau);

[ zdopt, Q, kara ] = BFGS(zd, tau, param, rho, opts);

zdopt
Q
kara

[ x, ~, ~ ] = solver(zdopt, param, var, rho);
wizualizacja(x, param);
