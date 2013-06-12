close all;
clear all;
format long e;

param;
%init = x0(LEO, 270, 0, 0, m);
init = x0(LEO, 230, 0, 3, m);
%init = [1 1 1 1 1];
u = [0; 0.2; 0; 0; pi; 0];
%u0 = [0 0; 0.2 pi; 0 0];
tau = [0 T/3 2*T/3 T];
ogr = [10 -10; 10 -10; 10 -10; -inf inf; -inf inf; -inf inf];
[ uopt, Q ] = BFGS(init, h, tau, u, ogr);

[ t, x, psi, grad ] = solver(init, h, tau, uopt);
wizualizacja(x);
%Q = solverSzybki(init, h, tau, u)