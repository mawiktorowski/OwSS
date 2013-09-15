function out = rhs(t, x, u, param)
% wyliczanie prawej strony rownania modelu dla RK4

xm = x(1) - param.D * cos(param.omega * t); 
ym = x(2) - param.D * sin(param.omega * t);
xm2 = xm * xm;
ym2 = ym * ym;
x12 = x(1) * x(1);
x22 = x(2) * x(2);
rE = realsqrt(x12 + x22); 
rM = realsqrt(xm2 + ym2);
rE3 = rE * rE * rE;
rM3 = rM * rM *rM;
coeff1 = param.gE / rE3;
coeff2 = param.gM / rM3;

out = zeros(1,5);
out(1) = x(3);
out(2) = x(4);
out(3) = - coeff1 * x(1) - coeff2 * xm + u(1) * cos(u(2)) / x(5);
out(4) = - coeff1 * x(2) - coeff2 * ym + u(1) * sin(u(2)) / x(5);
out(5) = param.C1 * u(1);
end
