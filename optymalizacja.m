close all;
clear all;
format long e;

param;
% kolejnosc zapisu w wektorze [theta alfa dV wszystkie u wszystkie psi]
%init = x0(LEO, 270, 0, 0, m);
%init = x0(LEO, 230, 0, 3, m); 
%init = [1 1 1 1 1];
thetaRad = 230 * 2 * pi / 360;
alphaRad = 0;
dVE = 3;
zd = [thetaRad; alphaRad; dVE; 0; 0.2; 0; 0; pi; 0];
u = [0; 0.2; 0; 0; pi; 0];
%u0 = [0 0; 0.2 pi; 0 0];
tau = [0 T/3 2*T/3 T];
ogr = [0 2*pi; 0 pi; 0 10; 0 10; 0 10; 0 10; -inf inf; -inf inf; -inf inf];

tic
[ zdopt, Q ] = BFGS(zd, h, tau, ogr);
toc

[ t, x, psi, grad ] = solver(zdopt, h, tau);
wizualizacja(x);
%Q = solverSzybki(init, h, tau, u)