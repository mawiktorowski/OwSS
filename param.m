global T h K1 K2
global mr
global D mu restmu
global aE aM LEO LMO rM VM
global C1 C2

% tymczasowe !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
global dV

% parametry symulacji
T = 1;
h = 0.0001;
mE = 5.9742e24;
mM = 7.3483e22;
mu = mM / (mM + mE);
restmu = 1 - mu;

% chwilowo tutaj, p�niej trzeba zrobi� jako parametr %%%%%%%%%%%%%%
%%% nieu�ywane !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dV = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K1 = 0; % warto�� wsp�czynnika kary dla organiczenia zwi�zanego z mas�
K2 = 1; % warto�� wsp�czynnika zewn�trznej funkcji kary

% parametry nieprzeskalowane

%G = 6.672e-20;
D = 384400;
%omega = sqrt(G * (mM + mE) / D^3);

aE = 6378;
aM = 1738;
LEO = 463;
LMO = 100;

rM = (LMO + aM) / D;
VM = sqrt(mu/rM);

%%% wsp�czynnik spalania
C1 = 0.1 / 300;
C2 = 0.1 / 3000;

% parametry modelu (rakieta)
mr = 10000 / 1000; % masa rakiety bez paliwa
m = 20000;