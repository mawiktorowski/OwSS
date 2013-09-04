close all;
clear all;
format long e;

[param, opts] = parametry();

thetaRad = 230 * 2 * pi / 360;
dVE = 3;
%u = [0; 0; 0; 0; 0; 0];
u = [0.2; 0.2; 0.2; 0; 0; 0];

zd = [thetaRad; dVE; u];

k = 10^6;
epKara = 10^-6;
c = 10;

T = 1;
tau = [0 T/3 2*T/3 T];
rho = 1;

var = obliczenia(param.h0, tau);

i = 1;
while rho < k
    [ zdopt, Q, kara ] = BFGS(zd, tau, param, rho, opts);
    disp([Q kara]);
    if (norm(zd - zdopt) < epKara)
       break 
    end
    zd = zdopt;
    rho = rho * c;
end

zdopt
Q
kara

[ x, ~, ~ ] = solver(zdopt, param, var, 0);
wizualizacja(x, param);