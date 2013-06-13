close all;
clear all;
format long;

param;
%init = x0(LEO, 270, 0, 0, m);
init = [1 1 1 1 10];
%init = [1 1 1 1 1];
ep=1e-6;
out = zeros(1,5);

Q0 = kosztSzybkiX0(init, [0 0; 1 pi], h, [0 T/2 T]);

for i=1:5
    x0tmp = init;
    x0tmp(i) = init(i)+ep;
    Q = kosztSzybkiX0(x0tmp, [0 0; 1 pi], h, [0 T/2 T]);
    out(i) = (Q - Q0) ./ ep;
end

[ t, x, psi] = solverX0(init, [0 0; 1 pi], h, [0 T/2 T]);
wizualizacja(x);

disp('Porównanie');
disp ([out', -psi(1,:)']);
disp('Procentowo (%)');
disp (abs((out' + psi(1,:)') ./ out') * 100);