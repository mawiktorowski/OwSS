function [param, BFGSopts] = parametry()
% funkcja ustawiajaca parametry

% dane liczbowe przed skalowaniem
aE = 6378;
aM = 1738;
D = 384400;
muE = 3.9860e5;
muM = 4.9028e3;
C1 = - 1/30;
C2 = - 1/15;
m0 = 2e5;
mr = 1e5;
dVmax = 6;
urmax = 1.5;

% obliczenia przed skalowaniem
mu = muM / (muE + muM);
restmu = 1 - mu;
omega = sqrt((muE + muM) / D^3);

% jednostki
DU = D;
TU = 1 / omega;
FU = urmax; 

% parametry ustawiane symulacji
h0 = 0.0001;
K1 = 1;
K2 = 0;
LEO = 463;
LMO = 100;

% parametry algorytmu BFGS
ep0 = 1e-6;     % minimalna norma gradientu
ep1 = 0.5e-12;  % odnowa algorytmu
ep2 = 1e-5;   	% odnowa algorytmu
epK0 = 1e-16; 	% dokladnosc kontrakcji - d == kierunek najszybszego spadku
epK1 = 1e-6;    % dokladnosc kontrakcji - d != kierunek najszybszego spadku
%MAX_ITER = 20000;  % ilosc iteracji
MAX_ITER = 20;  % ilosc iteracji
epOgr = 1e-10;

% parametry symulacji czesc obliczana

% dane liczbowe po przeskalowaniu
aE = aE/DU;
aM = aM/DU;
% D = D/DU;
% muE = muE * TU^2 / DU^3;
% muM = muM * TU^2 / DU^3;
% omega = omega * TU;
C1 = C1 * DU / TU;
C2 = C2 * DU / TU;
m0 = m0 * DU / FU / TU^2;
mr = mr * DU / FU / TU^2;
dVmax = dVmax * TU / DU;
% urmax = urmax / FU;

% wyliczane dane
rE = aE + (LEO / DU);
VE = sqrt(restmu/rE);
rM = aM + (LMO / DU);
VM = sqrt(mu/rM);
rM2 = rM * rM;
VM2 = VM * VM;
ogr = [0 2*pi; 0 dVmax; 0 1; -inf inf];

%na razie nie wiadomo co bedzie potrzebne
%param = struct('omega', omega, 'muE', muE, 'muM', muM, 'D', D, 'rE', rE, 'VE', VE, 'rM', rM, 'VM', VM, 'rM2', rM2, 'VM2', VM2, 'm0', m0, 'mr', mr, 'C1', C1, 'C2', C2, 'K1', K1, 'K2', K2, 'dVmax', dVmax);
param = struct('mu', mu, 'restmu', restmu, 'rE', rE, 'VE', VE, 'rM2', rM2, 'VM2', VM2, 'm0', m0, 'mr', mr, 'C1', C1, 'C2', C2, 'K1', K1, 'K2', K2, 'dVmax', dVmax, 'h0', h0, 'ogr', ogr);
BFGSopts = struct('ep0', ep0, 'ep1', ep1, 'ep2', ep2, 'epK0', epK0, 'epK1', epK1, 'MAX_ITER', MAX_ITER, 'epOgr', epOgr);
end
