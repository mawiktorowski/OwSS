% test Psi
close all;
clear all;
format long;

[param, ~] = parametry();

x0 = [1 1 1 1 1];
u = [0; 0; 0; 0];
h0 = param.h0;
T = 1;
tau = [0 T/2 T];
ep=1e-6;
rho = 1;

var = obliczenia(param.h0, tau);

[ dQ, H ] = sprawdzeniePsi(x0, u, param, var, ep, rho);

disp('Porównanie');
disp ([dQ, H]);
disp('Procentowo (%)');
disp (abs((dQ - H) ./ H) * 100);
