function out = rhsPsi(x, u, param)
% % wyliczanie prawej strony rownan w postaci kanonicznej

xm = x(3) - x(1); 
ym = x(4) - x(2);
xm2 = xm * xm;
ym2 = ym * ym;
x32 = x(3) * x(3);
x42 = x(4) * x(4);
rE = realsqrt(x32 + x42); 
rM = realsqrt(xm2 + ym2);
rE2 = rE * rE;
rM2 = rM * rM;
rE3 = rE2 * rE;
rM3 = rM2 * rM;
coeff1 = param.gE / param.D3;
coeff2 = param.gE / rE3;
coeff3 = param.gM / rM3;
coeff4 = coeff2 / rE2;
coeff5 = coeff3 / rM2;
x92 = x(9) * x(9);
sphi = sin(u(2));
cphi = cos(u(2));

d71 = - coeff5 * (rM2 - 3 * xm2);
d81 = coeff5 * 3 * xm * ym;
d72 = d81;
d82 = - coeff5 * (rM2 - 3 * ym2);
d73 = coeff4 * (rE2 - 3 * x32) + coeff5 * (rM2 - 3 * xm2);
d83 = - 3 * (coeff4 * x(3) * x(4) + coeff5 * xm * ym);
d74 = d83;
d84 = coeff4 * (rE2 - 3 * x42) + coeff5 * (rM2 - 3 * ym2);

out = zeros(1,20);
out(1)  = x(5);
out(2)  = x(6);
out(3)  = x(7);
out(4)  = x(8);
out(5)  = - coeff1 * x(1);
out(6)  = - coeff1 * x(2);
out(7)  = - coeff2 * x(3) - coeff3 * xm + u(1) * cphi / x(9);
out(8)  = - coeff2 * x(4) - coeff3 * ym + u(1) * sphi / x(9);
out(9)  = param.C1 * u(1);
% rownania sprzezone
out(10) = coeff1 * x(14) + d71 * x(16) + d81 * x(17);
out(11) = coeff1 * x(15) + d72 * x(16) + d82 * x(17);
out(12) = d73 * x(16) + d83 * x(17);
out(13) = d74 * x(16) + d84 * x(17);
out(14) = - x(10);
out(15) = - x(11);
out(16) = - x(12);
out(17) = - x(13);
out(18) = u(1) * (x(16) * cphi + x(17) * sphi) / x92;
% dH/du
out(19) = (x(16) * cphi + x(17) * sphi) / x(9) + x(18) * param.C1;
out(20) = u(1) * (- x(16) * sphi + x(17) * cphi) / x(9);
end
