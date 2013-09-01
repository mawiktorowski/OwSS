function out = rhs(x, u, param)
% wyliczanie prawej strony rownania modelu dla RK4

x1mu = x(1) + param.mu;
x1rmu = x(1) - param.restmu;
x22 = x(2) * x(2);
x1mu2 = x1mu * x1mu;
x1rmu2 = x1rmu * x1rmu;
rE = realsqrt(x1mu2 + x22); 
rM = realsqrt(x1rmu2 + x22);
rE3 = rE * rE * rE;
rM3 = rM * rM *rM;
coeff1 = param.restmu / rE3;
coeff2 = param.mu / rM3;

out = zeros(1,5);
out(1) = x(3);
out(2) = x(4);
out(3) = x(1) + 2 * x(4) - coeff1 * x1mu - coeff2 * x1rmu + u(1) * cos(u(2)) / x(5);
out(4) = x(2) - 2 * x(3) - coeff1 * x(2) - coeff2 * x(2)  + u(1) * sin(u(2)) / x(5);
out(5) = param.C1 * u(1);
end
