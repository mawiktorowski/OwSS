global T h K1 K2
global mr m
global D mu restmu
global aE aM LEO LMO rM VM rE VE
global C1 C2
global ep0 ep1 ep2 epK0 epK1 MAX_ITER

% tymczasowe !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
global dV

% parametry symulacji
T = 1;
h = 0.0001;
mE = 5.9742e24;
mM = 7.3483e22;
mu = mM / (mM + mE);
restmu = 1 - mu;

% chwilowo tutaj, póŸniej trzeba zrobiæ jako parametr %%%%%%%%%%%%%%
%%% nieu¿ywane !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dV = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K1 = 0; % wartoœæ wspó³czynnika kary dla organiczenia zwi¹zanego z mas¹
K2 = 1; % wartoœæ wspó³czynnika zewnêtrznej funkcji kary

% parametry nieprzeskalowane

%G = 6.672e-20;
D = 384400;
%omega = sqrt(G * (mM + mE) / D^3);

aE = 6378;
aM = 1738;
LEO = 463;
LMO = 100;

rE = (aE + LEO) / D;
VE = sqrt(restmu/rE);

rM = (LMO + aM) / D;
VM = sqrt(mu/rM);

%%% wspó³czynnik spalania
C1 = 0.1 / 300;
C2 = 0.1 / 3000;

% parametry modelu (rakieta)
mr = 10000 / 1000; % masa rakiety bez paliwa
m = 20000 / 1000;

% parametry algorytmu BFGS
ep0 = 1e-8;     % minimalna norma gradientu
ep1 = 0.5e-16;  % odnowa algorytmu
ep2 = 0.5;   	% odnowa algorytmu
epK0 = 1e-16; 	% dokladnosc kontrakcji - d == kierunek najszybszego spadku
epK1 = 1e-6;    % dokladnosc kontrakcji - d != kierunek najszybszego spadku
MAX_ITER = 10;

%opts = struct('ep0', 1e-8, 'ep1', 0.5e-16, 'ep2', 0.5, 'epK0', 1e-16, 'epK1', 1e-6, 'MAX_ITER', 10);





