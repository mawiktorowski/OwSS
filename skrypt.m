close all;
clear all;

%format long e;

param;
thetaRad = 270 * 2 * pi / 360;
alphaRad = 0;
dVE = 0;

%init = x0(LEO, 270, 0, 0, m);
%init = x0(LEO, 290, 0, 3, m);
%init = [1 1 1 1 20];
zd = [thetaRad; alphaRad; dVE; 0; 0];

%init = x0(LEO, 270, 0, 0, m)
%init = [ D+rM 0 0 VM + D * omegaM m];
%init = x0(LEO, 270, 0, 0, m);
%init = x0(LEO, -116, 0, 3.2, m);
% % solver rk4
%[ t, x, psi ] = solver(init, h, [0 T/2 T], [0 0; 1 pi]);
[ t, x, psi ] = solver(zd, h, [0 T]);
wizualizacja(x);
% %y = ode4(@(t,x) rhs(x, u(t))',t,init);
% 
% %delta = x - y;
% 
% % solver rk4 równañ sprzê¿onych w ty³
% xT = x(end,:);
% initBack = canT(xT, t(end));
% [tCan,xCan] = rk4Back(@rhsCan, @u, h, T, initBack);
% 
% %yCan = ode4(@(t,x) rhsCan(x, u(t))',fliplr(tCan),xT);
% 
% %yCan = flipud(yCan);
% 
% %deltaCan = xCan - yCan;
% 
% % sprawdzanie poprawnoœci równañ