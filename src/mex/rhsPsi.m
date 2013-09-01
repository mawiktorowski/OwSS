function out = rhsPsi(x, u, param)
% % wyliczanie prawej strony rownan w postaci kanonicznej

x1mu = x(1) + param.mu;
x1rmu = x(1) - param.restmu;
x22 = x(2) * x(2);
x1mu2 = x1mu * x1mu;
x1rmu2 = x1rmu * x1rmu;
rE = realsqrt(x1mu2 + x22); 
rM = realsqrt(x1rmu2 + x22);
rE2 = rE * rE;
rM2 = rM * rM;
rE3 = rE2 * rE;
rM3 = rM2 * rM;
coeff1 = param.restmu / rE3;
coeff2 = param.mu / rM3;
coeff3 = coeff1 / rE2;
coeff4 = coeff2 / rM2;

x52 = x(5) * x(5);
sphi = sin(u(2));
cphi = cos(u(2));

d31 = -1 + coeff3 * (rE2 - 3 * x1mu2) + coeff4 * (rM2 - 3 * x1rmu2);
d41 = - 3 * coeff3 * x1mu * x(2) - 3 * coeff4 * x1rmu * x(2);
d32 = d41;
d42 = -1 + coeff3 * (rE2 - 3 * x22) + coeff4 * (rM2 - 3 * x22);

out = zeros(1,12);
out(1) = x(3);
out(2) = x(4);
out(3) = x(1) + 2 * x(4) - coeff1 * x1mu - coeff2 * x1rmu + u(1) * cphi / x(5);
out(4) = x(2) - 2 * x(3) - coeff1 * x(2) - coeff2 * x(2)  + u(1) * sphi / x(5);
out(5) = param.C1 * u(1);
% rownania sprzezone
out(6) = x(8) * d31 + x(9) * d41;
out(7) = x(8) * d32 + x(9) * d42;
out(8) = - x(6) + 2 * x(9);
out(9) = - x(7) - 2 * x(8);
out(10) = x(8) * u(1) * cphi / x52 + x(9) * u(1) * sphi / x52;
% dH/du
out(11) = x(8) * cphi / x(5) + x(9) * sphi / x(5) + x(10) * param.C1;
out(12) = - x(8) * u(1) * sphi / x(5) + x(9) * u(1) * cphi / x(5);
end
