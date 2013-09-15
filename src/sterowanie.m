function sterowanie(zd, tau, param)
% rysowanie wykresów sterowañ

lzd = length(zd);
lu = (lzd - 2)/2;
urend = lu + 2;
phistart = urend + 1;

theta = radtodeg(wrapTo2Pi(zd(1)));
v0 = zd(2) * param.DU / param.TU;
ur = zd(3:urend) * param.FU;
phi = radtodeg(wrapToPi(zd(phistart:end)));

fprintf(['K¹t startowy: ', num2str(theta) ' (deg)\n']);
fprintf(['Wartoœæ impulsu startowego: ', num2str(v0) ' (km/s)\n']);

tauSec = tau * param.TU;

figure;
hold on;
stairs(tauSec, [ur; 0]);
hold off;
xlabel('T(s)');
ylabel('ur (kN)');
title('Sterowanie ci¹gniem w funkcji czasu');

figure;
hold on;
stairs(tauSec, [phi; 0]);
hold off;
xlabel('T(s)');
ylabel('phi (deg)');
title('Sterowanie k¹tem odchylenia dyszy w funkcji czasu');

end
