function [ Q, kara ] = kosztSzybki(zd, param, var, rho)
% szybkie wyliczanie wskaznika jakosci

x = zeros(1,9);
us = zeros(2,1);
u = zd(3:end);

x(1) = param.D;
x(2) = 0;
x(3) = param.rE * cos(zd(1));
x(4) = param.rE * sin(zd(1));
x(5) = 0;
x(6) = param.omega * param.D;
x(7) = - (param.VE + zd(2)) * sin(zd(1));
x(8) = (param.VE + zd(2)) * cos(zd(1));
x(9) = param.m0 * exp(param.C2 * zd(2));

for j = 1:var.ldtau
    jj = var.ldtau + j;
    us(1) = u(j);
    us(2) = u(jj);
    hj = var.h(j);
    h2j = var.h2(j);
    h3j = var.h3(j);
    h6j = var.h6(j);
    for i = var.cn(j):var.cn(j+1)-1
        dx1 = rhs(x, us, param);
        farg = x + h2j * dx1;
        dx2 = rhs(farg, us, param);
        farg = x + h2j * dx2;
        dx3 = rhs(farg, us, param);
        farg = x + hj * dx3;
        dx4 = rhs(farg, us, param);
        x = x + h3j * (dx2 + dx3) + h6j * (dx1 + dx4);
    end
end

xm = x(3) - x(1);
ym = x(4) - x(2);
um = x(7) - x(5);
vm = x(8) - x(6);

xm2 = xm * xm;
ym2 = ym * ym;
um2 = um * um;
vm2 = vm * vm;

beta1 = xm2 + ym2 - param.rM2;
beta2 = um2 + vm2 - param.VM2;
beta3 = xm * um + ym * vm;

k1 = 0.25 * beta1 * beta1;
k2 = 0.25 * beta2 * beta2;
k3 = 0.5 * beta3 * beta3;
if x(9) < param.mr
    k4 = 0.5 * (x(9) - param.mr)^2;
else
    k4 = 0;
end

kara = k1 + k2 + k3 + k4;
Q = - param.K1 * x(9) + param.K2 * var.tf + rho * kara;

end