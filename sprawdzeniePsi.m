close all;
clear all;
format long;

param;
%init = x0(LEO, 270, 0, 0, m);
init = x0(LEO, 230, 0, 3, m);
%init = [1 1 1 1 1];
ep=1e-6;
out = zeros(1,5);

Q0 = solverSzybki(init, h, [0 T/2 T], [0 0; 1 pi]);

for i=1:5
    x0tmp = init;
    x0tmp(i) = init(i)+ep;
    Q = solverSzybki(x0tmp, h, [0 T/2 T], [0 0; 1 pi]);
    out(i) = (Q - Q0) ./ ep;
end

[ t, x, psi, grad ] = solver(init, h, [0 T/2 T], [0 0; 1 pi]);
wizualizacja(x);

disp('Porównanie');
disp ([out', -psi(1,:)']);
disp('Procentowo (%)');
disp (abs((out' + psi(1,:)') ./ out') * 100);