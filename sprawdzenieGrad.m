close all;
clear all;
format long e;

param;
%init = x0(LEO, 270, 0, 0, m);
init = x0(LEO, 290, 0, 3, m);
%init = [1 1 1 1 20];
u0 = [0 0.2 0 0 pi 0];
%u0 = [0 0; 0.2 pi; 0 0];
tau = [0 T/3 2*T/3 T];
ep=1e-6;
out = zeros((2 * (length(tau) - 1)),1);

Q0 = solverSzybki(init, h, tau, u0);

for i=1:(2 * (length(tau) - 1))
    u0tmp = u0;
    u0tmp(i) = u0(i)+ep;
    Q = solverSzybki(init, h, [0 T/3 2*T/3 T], u0tmp);
    out(i) = (Q - Q0) ./ ep;
end

[ ~, grad ] = solverSzybki(init, h, [0 T/3 2*T/3 T], u0);

disp('Porównanie');
disp ([out grad]);
disp('Procentowo (%)');
disp (abs([out - grad] ./ out) * 100);