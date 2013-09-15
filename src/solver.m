function [ x, Q, psi, grad ] = solver(zd, param, var, rho)
% solver modelu rozwiazujacy rownania w przod i wstecz z zapamietywaniem

x = zeros(var.iter,9);
us = zeros(2,1);
u = zd(3:end);

x(1,1) = param.D;
x(1,2) = 0;
x(1,3) = param.rE * cos(zd(1));
x(1,4) = param.rE * sin(zd(1));
x(1,5) = 0;
x(1,6) = param.omega * param.D;
x(1,7) = - (param.VE + zd(2)) * sin(zd(1));
x(1,8) = (param.VE + zd(2)) * cos(zd(1));
x(1,9) = param.m0 * exp(param.C2 * zd(2));

for j = 1:var.ldtau
    jj = var.ldtau + j;
    us(1) = u(j);
    us(2) = u(jj);
    hj = var.h(j);
    h2j = var.h2(j);
    h3j = var.h3(j);
    h6j = var.h6(j);
    for i = var.cn(j):var.cn(j+1)-1
        xc = x(i,:);
        dx1 = rhs(xc, us, param);
        farg = xc + h2j * dx1;
        dx2 = rhs(farg, us, param);
        farg = xc + h2j * dx2;
        dx3 = rhs(farg, us, param);
        farg = xc + hj * dx3;
        dx4 = rhs(farg, us, param);
        x(i+1,:) = xc + h3j * (dx2 + dx3) + h6j * (dx1 + dx4);
    end
end

xT = x(var.iter,:);

xm = xT(3) - xT(1);
ym = xT(4) - xT(2);
um = xT(7) - xT(5);
vm = xT(8) - xT(6);

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

if x(9) <= param.mr
    diffK4 = x(9) - param.mr;
else
    diffK4 = 0;
end

psiT = zeros(1,9);

psiT(1) =   rho * ( beta1 * xm + beta3 * um);
psiT(2) =   rho * ( beta1 * ym + beta3 * vm);
psiT(3) = - psiT(1);
psiT(4) = - psiT(2);
psiT(5) =   rho * ( beta2 * um + beta3 * xm);
psiT(6) =   rho * ( beta2 * vm + beta3 * ym);
psiT(7) = - psiT(5);
psiT(8) = - psiT(6);
psiT(9) = param.K1 - rho * diffK4;

x = [x zeros(var.iter, 7)];

x(var.iter, 10:20) = [psiT 0 0];

for j = var.ldtau:-1:1
    jj = var.ldtau + j;
    us(1) = u(j);
    us(2) = u(jj);
    hj = var.h(j);
    h2j = var.h2(j);
    h3j = var.h3(j);
    h6j = var.h6(j);
    for i = (var.cn(j+1)-1):-1:var.cn(j)
        xc = x(i+1,:);
        dx1 = rhsPsi(xc, us, param);
        farg = xc - h2j * dx1;
        dx2 = rhsPsi(farg, us, param);
        farg = xc - h2j * dx2;
        dx3 = rhsPsi(farg, us, param);
        farg = xc - hj * dx3;
        dx4 = rhsPsi(farg, us, param);
        x(i,:) = xc - h3j * (dx2 + dx3) - h6j * (dx1 + dx4);
    end
    gradU(j) = x(i,11);
    gradU(jj) = x(i,12);
    x(i, 11:12) = [0 0]; % gamma(t_i+1) = 0
end

dx3theta = - param.rE * sin(zd(1));
dx4theta = param.rE * cos(zd(1));
dx7theta = - (zd(2) + param.VE) * cos(zd(1));
dx8theta = - (zd(2) + param.VE) * sin(zd(1));
dx7dV = - sin(zd(1));
dx8dV = cos(zd(1));
dx9dV = param.m0 * param.C2 * exp(param.C2 * zd(2));

grad = zeros(length(zd),1);
grad(1) = - x(1,12) * dx3theta - x(1,13) * dx4theta - x(1,16) * dx7theta - x(1,17) * dx8theta;
grad(2) = - x(1,16) * dx7dV - x(1,17) * dx8dV - x(1,18) * dx9dV;
grad(3:end) = gradU;

psi = x(:,10:20);
x = x(:,1:9);
end
