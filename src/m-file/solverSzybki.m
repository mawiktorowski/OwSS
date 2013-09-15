function [ Q, grad ] = solverSzybki(zd, param, var, rho)
% solver szybki rozwiazujacy rownania zwraca wskaznik jakosci i gradient

x = zeros(1,5);
gradU = zeros(2 * var.ldtau,1);
us = zeros(2,1);
u = zd(3:end);
t = 0;

x(1) = param.rE * cos(zd(1));
x(2) = param.rE * sin(zd(1));
x(3) = - (param.VE + zd(2)) * sin(zd(1));
x(4) = (param.VE + zd(2)) * cos(zd(1));
x(5) = param.m0 * exp(param.C2 * zd(2));

for j = 1:var.ldtau
    jj = var.ldtau + j;
    us(1) = u(j);
    us(2) = u(jj);
    hj = var.h(j);
    h2j = var.h2(j);
    h3j = var.h3(j);
    h6j = var.h6(j);
    for i = var.cn(j):var.cn(j+1)-1
        dx1 = rhs(t, x, us, param);
        targ = t + h2j;
        farg = x + h2j * dx1;
        dx2 = rhs(targ, farg, us, param);
        farg = x + h2j * dx2;
        dx3 = rhs(targ, farg, us, param);
        targ = t + hj;
        farg = x + hj * dx3;
        dx4 = rhs(targ, farg, us, param);
        x = x + h3j * (dx2 + dx3) + h6j * (dx1 + dx4);
        t = t + hj;
    end
end

xm = x(1) -               param.D * cos(param.omega * t);
ym = x(2) -               param.D * sin(param.omega * t);
um = x(3) + param.omega * param.D * sin(param.omega * t);
vm = x(4) - param.omega * param.D * cos(param.omega * t);

xm2 = xm * xm;
ym2 = ym * ym;
um2 = um * um;
vm2 = vm * vm;

beta1 = xm2 + ym2 - param.rM2;
beta2 = um2 + vm2 - param.VM2;
beta3 = xm * um + ym * vm;

if x(5) <= param.mr
    diffK4 = x(5) - param.mr;
else
    diffK4 = 0;
end

psi = zeros(1,5);

psi(1) =  - rho * ( beta1 * xm + beta3 * um);
psi(2) =  - rho * ( beta1 * ym + beta3 * vm);
psi(3) =  - rho * ( beta2 * um + beta3 * xm);
psi(4) =  - rho * ( beta2 * vm + beta3 * ym);
psi(5) = param.K1 - rho * diffK4;

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
    gradU(j) = x(11);
    gradU(jj) = x(12);
    x(11:12) = [0 0]; % gamma(t_i+1) = 0
end

dx1theta = - param.rE * sin(zd(1));
dx2theta = param.rE * cos(zd(1));
dx3theta = - (zd(2) + param.VE) * cos(zd(1));
dx4theta = - (zd(2) + param.VE) * sin(zd(1));
dx3dV = - sin(zd(1));
dx4dV = cos(zd(1));
dx5dV = param.m0 * param.C2 * exp(param.C2 * zd(2));

grad = zeros(length(zd),1);
grad(1) = - x(6) * dx1theta - x(7) * dx2theta - x(8) * dx3theta - x(9) * dx4theta;
grad(2) = - x(8) * dx3dV - x(9) * dx4dV - x(10) * dx5dV;
grad(3:end) = gradU;

end
