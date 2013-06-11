close all;
clear all;

param;
init = x0(LEO, 270, 0, 2, m)
%init = [ D+rM 0 0 VM + D * omegaM m];
% % solver rk4
%[t,x] = rk4(@rhs, @u, h, T, init);
[ t, x, psi ] = solver(init, h, [0 T/2 T], [0 0; 1 pi]);
wizualizacja(x,t);
wizualizacjaPrim(x,t);
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