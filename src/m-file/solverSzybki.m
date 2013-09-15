function [ Q, grad ] = solverSzybki(zd, param, var, rho)
% solver szybki rozwiazujacy rownania zwraca wskaznik jakosci i gradient

x = zeros(1,9);
gradU = zeros(2 * var.ldtau,1);
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

if x(9) <= param.mr
    diffK4 = x(9) - param.mr;
else
    diffK4 = 0;
end

psi = zeros(1,9);

psi(1) =   rho * ( beta1 * xm + beta3 * um);
psi(2) =   rho * ( beta1 * ym + beta3 * vm);
psi(3) = - psi(1) ;
psi(4) = - psi(2);
psi(5) =   rho * ( beta2 * um + beta3 * xm);
psi(6) =   rho * ( beta2 * vm + beta3 * ym);
psi(7) = - psi(5);
psi(8) = - psi(6);
psi(9) = param.K1 - rho * diffK4;

x = [x psi 0 0];

for j = var.ldtau:-1:1
    jj = var.ldtau + j;
    us(1) = u(j);
    us(2) = u(jj);
    hj = var.h(j);
    h2j = var.h2(j);
    h3j = var.h3(j);
    h6j = var.h6(j);
    for i = (var.cn(j+1)-1):-1:var.cn(j)
        dx1 = rhsPsi(x, us, param);
        farg = x - h2j * dx1;
        dx2 = rhsPsi(farg, us, param);
        farg = x - h2j * dx2;
        dx3 = rhsPsi(farg, us, param);
        farg = x - hj * dx3;
        dx4 = rhsPsi(farg, us, param);
        x = x - h3j * (dx2 + dx3) - h6j * (dx1 + dx4);
    end
    gradU(j) = x(19);
    gradU(jj) = x(20);
    x(19:20) = [0 0]; % gamma(t_i+1) = 0
end

dx3theta = - param.rE * sin(zd(1));
dx4theta = param.rE * cos(zd(1));
dx7theta = - (zd(2) + param.VE) * cos(zd(1));
dx8theta = - (zd(2) + param.VE) * sin(zd(1));
dx7dV = - sin(zd(1));
dx8dV = cos(zd(1));
dx9dV = param.m0 * param.C2 * exp(param.C2 * zd(2));

grad = zeros(length(zd),1);
grad(1) = - x(12) * dx3theta - x(13) * dx4theta - x(16) * dx7theta - x(17) * dx8theta;
grad(2) = - x(16) * dx7dV - x(17) * dx8dV - x(18) * dx9dV;
grad(3:end) = gradU;

end
