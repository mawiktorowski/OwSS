function [param, opts] = parametry()
% funkcja ustawiajaca parametry

% dane liczbowe przed skalowaniem
aE = 6378;
aM = 1738;
D = 384400;
gE = 3.9860e5;
gM = 4.9028e3;
C1 = - 1/30;
C2 = - 1/15;
m0 = 2e5;
mr = 1e5;
dVmax = 6;
urmax = 1.5;

% obliczenia przed skalowaniem
omega = sqrt(gE / D^3);

% jednostki
DU = aE;
TU = sqrt(aE^3 / gE);
FU = urmax; 

% parametry ustawiane symulacji
%h0 = 0.0001;
h0 = 0.01;
K1 = 1;
K2 = 1;
LEO = 463;
LMO = 100;

% parametry algorytmu BFGS
ep0 = 1e-8;     % minimalna norma gradientu
ep1 = 0.5e-16;  % odnowa algorytmu
ep2 = 1e-7;   	% odnowa algorytmu
epK0 = 1e-16; 	% dokladnosc kontrakcji - d == kierunek najszybszego spadku
epK1 = 1e-6;    % dokladnosc kontrakcji - d != kierunek najszybszego spadku
%MAX_ITER = 20000;  % ilosc iteracji
MAX_ITER = 1000;  % ilosc iteracji
epOgr = 1e-10;

%parametry algorytmu funkcji kary
epKara = 1e-1;
k = 125; %ograniczenie wspó³czynnika kary
c = 2; %mno¿nik wspó³czynnika kary

% parametry symulacji czesc obliczana

% dane liczbowe po przeskalowaniu
aE = aE/DU;
aM = aM/DU;
D = D/DU;
gE = gE * TU^2 / DU^3;
gM = gM * TU^2 / DU^3;
omega = omega * TU;
C1 = C1 * DU / TU;
C2 = C2 * DU / TU;
m0 = m0 * DU / FU / TU^2;
mr = mr * DU / FU / TU^2;
dVmax = dVmax * TU / DU;
% urmax = urmax / FU;

% wyliczane dane
D3 = D * D * D;
rE = aE + (LEO / DU);
VE = sqrt(gE / rE);
rM = aM + (LMO / DU);
VM = sqrt(gM / rM);
rM2 = rM * rM;
VM2 = VM * VM;
ogr = [0 2*pi; 0 dVmax; 0 1; -inf inf];

%na razie nie wiadomo co bedzie potrzebne
%param = struct('omega', omega, 'muE', muE, 'muM', muM, 'D', D, 'rE', rE, 'VE', VE, 'rM', rM, 'VM', VM, 'rM2', rM2, 'VM2', VM2, 'm0', m0, 'mr', mr, 'C1', C1, 'C2', C2, 'K1', K1, 'K2', K2, 'dVmax', dVmax);
param = struct('omega', omega, 'gE', gE, 'gM', gM, 'D', D, 'D3', D3, 'rE', rE, 'VE', VE, 'rM2', rM2, 'VM2', VM2, 'm0', m0, 'mr', mr, 'C1', C1, 'C2', C2, 'K1', K1, 'K2', K2, 'dVmax', dVmax, 'h0', h0, 'ogr', ogr, 'DU', DU, 'TU', TU, 'FU', FU);
opts = struct('ep0', ep0, 'ep1', ep1, 'ep2', ep2, 'epK0', epK0, 'epK1', epK1, 'MAX_ITER', MAX_ITER, 'epOgr', epOgr, 'epKara', epKara, 'c', c, 'k', k);
end
