close all;
clear all;
format long;

param;
%init = x0(LEO, 270, 0, 0, m);
init = x0(LEO, 290, 0, 3, m);
u0 = [0 0; 0.2 pi; 0 0];
tau = [0 T/3 2*T/3 T];
ep=1e-6;
out = zeros(length(tau)-1,2);

Q0 = kosztSzybki(init, h, tau, u0);

for j=1:2
    for i=1:length(tau)-1
        u0tmp = u0;
        u0tmp(i,j) = u0(i,j)+ep;
        Q = kosztSzybki(init, h, [0 T/3 2*T/3 T], u0tmp);
        out(i,j) = (Q - Q0) ./ ep;
    end
end

[ t, x, psi, grad ] = solver(init, h, [0 T/3 2*T/3 T], u0);
wizualizacja(x,t);

disp('Porównanie');
disp (out);
disp (grad);
disp('Procentowo');
disp ([out - grad] ./ out);