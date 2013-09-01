% test Psi
close all;
clear all;
format long;

[param, ~] = parametry();

x0 = [1 1 1 1 1];
psi0 = [1 1 1 1 1];
u = [1;1];

ep=1e-6;

[ dH, psi ] = sprawdzenieH(x0, psi0, u, param, ep);

disp('Porównanie');
disp ([dH, psi]);
disp('Procentowo (%)');
disp (abs((dH - psi) ./ psi) * 100);
