function out = rhsPsi(t, x, u, param)
% % wyliczanie prawej strony rownan w postaci kanonicznej

xm = x(1) - param.D * cos(param.omega * t); 
ym = x(2) - param.D * sin(param.omega * t);
xm2 = xm * xm;
ym2 = ym * ym;
x12 = x(1) * x(1);
x22 = x(2) * x(2);
rE = realsqrt(x12 + x22); 
rM = realsqrt(xm2 + ym2);
rE2 = rE * rE;
rM2 = rM * rM;
rE3 = rE2 * rE;
rM3 = rM2 * rM;
coeff1 = param.gE / rE3;
coeff2 = param.gM / rM3;
coeff3 = coeff1 / rE2;
coeff4 = coeff2 / rM2;
x52 = x(5) * x(5);
sphi = sin(u(2));
cphi = cos(u(2));

d31 = coeff3 * (rE2 - 3 * x12) + coeff4 * (rM2 - 3 * xm2);
d41 = - 3 * (coeff3 * x(1) * x(2) + coeff4 * xm * ym);
d32 = d41;
d42 = coeff3 * (rE2 - 3 * x22) + coeff4 * (rM2 - 3 * ym2);

out = zeros(1,12);
out(1) = x(3);
out(2) = x(4);
out(3) = - coeff1 * x(1) - coeff2 * xm + u(1) * cos(u(2)) / x(5);
out(4) = - coeff1 * x(2) - coeff2 * ym + u(1) * sin(u(2)) / x(5);
out(5) = param.C1 * u(1);
% rownania sprzezone
out(6) = d31 * x(8) + d41 * x(9);
out(7) = d32 * x(8) + d42 * x(9);
out(8) = - x(6);
out(9) = - x(7);
out(10) = u(1) * (x(8) * cphi + x(9) * sphi) / x52;
% dH/du
out(11) = (x(8) * cphi + x(9) * sphi) / x(5) + x(10) * param.C1;
out(12) = u(1) * (- x(8) * sphi + x(9) * cphi) / x(5);
end
