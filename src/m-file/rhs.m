function out = rhs(x, u, param)
% wyliczanie prawej strony rownania modelu dla RK4

xm = x(3) - x(1); 
ym = x(4) - x(2);
xm2 = xm * xm;
ym2 = ym * ym;
x32 = x(3) * x(3);
x42 = x(4) * x(4);
rE = realsqrt(x32 + x42); 
rM = realsqrt(xm2 + ym2);
rE3 = rE * rE * rE;
rM3 = rM * rM *rM;
coeff1 = param.gE / param.D3;
coeff2 = param.gE / rE3;
coeff3 = param.gM / rM3;

out = zeros(1,9);
out(1) = x(5);
out(2) = x(6);
out(3) = x(7);
out(4) = x(8);
out(5) = - coeff1 * x(1);
out(6) = - coeff1 * x(2);
out(7) = - coeff2 * x(3) - coeff3 * xm + u(1) * cos(u(2)) / x(9);
out(8) = - coeff2 * x(4) - coeff3 * ym + u(1) * sin(u(2)) / x(9);
out(9) = param.C1 * u(1);
end
