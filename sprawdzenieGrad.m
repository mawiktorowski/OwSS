close all;
clear all;
format long e;

param;
thetaRad = 290 * 2 * pi / 360;
alphaRad = 0;
dVE = 3;

%init = x0(LEO, 270, 0, 0, m);
%init = x0(LEO, 290, 0, 3, m);
%init = [1 1 1 1 20];
zd = [thetaRad; alphaRad; dVE; 0; 0.2; 0; 0; pi; 0];
%u0 = [0 0; 0.2 pi; 0 0];
tau = [0 T/3 2*T/3 T];
ep=1e-6;
out = zeros((2 * (length(tau) - 1)),1);

Q0 = solverSzybki(zd, h, tau);

for i=1:length(zd)
    zdtmp = zd;
    zdtmp(i) = zd(i)+ep;
    Q = solverSzybki(zdtmp, h, tau);
    out(i) = (Q - Q0) ./ ep;
end

[ ~, grad ] = solverSzybki(zd, h, tau);

[ t, x, psi, grad1 ] = solver(zd, h, tau);
wizualizacja(x);


disp('Porównanie');
disp ([out grad]);
disp('Procentowo (%)');
disp (abs([out - grad] ./ out) * 100);