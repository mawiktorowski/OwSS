function [ Q, grad ] = solverSzybki(zd, h0, tau)
% solver szybki modelu rozwiazujacy rownania bez zapamietywania trajektorii
% zwraca wskaznik jakosci i gradient

global mu restmu K2 rM VM mr
global C2 rE VE m

dtau = diff(tau);
n = ceil(dtau./h0);
h = dtau./n;
h2 = h/2;
h3 = h/3;
h6 = h/6;
cn = cumsum([1 n]);

x = zeros(1,5);
gradU = zeros(2*length(dtau),1);
us = zeros(2,1);
u = zd(4:end);

x(1) = rE * cos(zd(1)) - mu;
x(2) = rE * sin(zd(1));
x(3) = zd(3) * sin(zd(2) - zd(1)) - VE * sin(zd(1)) + x(2);
x(4) = zd(3) * cos(zd(2) - zd(1)) + VE * cos(zd(1)) - x(1);
x(5) = m * exp(- C2 * zd(3));
t = 0;

for j = 1:length(dtau)
    us(1) = u(j);
    us(2) = u(length(dtau) + j);
    f = @(x,u) rhs(x,u);
    for i = cn(j):cn(j+1)-1
        dx1 = f(x, us);
        dx2 = f(x + h2(j) * dx1, us);
        dx3 = f(x + h2(j) * dx2, us);
        dx4 = f(x + h(j) * dx3, us);
        x = x + h3(j) * (dx2 + dx3) + h6(j) * (dx1 + dx4);
        t = t + h(j);
    end
end

Q = koszt(x, t);

beta1 = rM^2 - (x(1) - restmu)^2 - x(2)^2;
beta2 = VM^2 - (x(3) - x(2))^2 - (x(4) + x(1))^2;
beta3 = - (x(1) - restmu) * (x(3) - x(2)) - x(2) * (x(4) + x(1));

    function y = diffK4
        if x(5) <= mr
            y = x(5) - mr;
        else
            y = 0;
        end
    end

psiT(1) = K2 * (beta1 * (x(1) - restmu) + beta2 * (x(4) + x(1)) + beta3 * x(3));
psiT(2) = K2 * (beta1 * x(2)            + beta2 * (x(2) - x(3)) + beta3 * (restmu + x(4)));
psiT(3) = K2 * (                          beta2 * (x(3) - x(2)) + beta3 * (x(1) - restmu));
psiT(4) = K2 * (                          beta2 * (x(4) + x(1)) + beta3 * x(2));
psiT(5) = 1 - K2 * diffK4;

x = [x psiT 0 0];

for j = length(dtau):-1:1
    us(1) = u(j);
    us(2) = u(length(dtau) + j);
    g = @(x,u) rhsGrad(x,u);
    for i = (cn(j+1)-1):-1:cn(j)
        dx1 = g(x, us);
        dx2 = g(x - h2(j) * dx1, us);
        dx3 = g(x - h2(j) * dx2, us);
        dx4 = g(x - h(j) * dx3, us);
        x = x - h3(j) * (dx2 + dx3) - h6(j) * (dx1 + dx4);
    end
    gradU(j) = x(11);
    gradU(length(dtau)+j) = x(12);
    x(11:12) = [0 0]; % gamma(t_i+1) = 0
end

dx1theta = - rE * sin(zd(1));
dx2theta = rE * cos(zd(1));
dx3theta = - zd(3) * cos(zd(2) - zd(1)) - VE * cos(zd(1)) + rE * cos(zd(1));
dx4theta = zd(3) * sin(zd(2) - zd(1)) - VE * sin(zd(1)) + rE * sin(zd(1));
dx3alpha = zd(3) * cos(zd(2) - zd(1));
dx4alpha = - zd(3) * sin(zd(2) - zd(1));
dx3dV = sin(zd(2) - zd(1));
dx4dV = cos(zd(2) - zd(1));
dx5dV = - m * C2 * exp(- C2 * zd(3));

grad = zeros(length(zd),1);
grad(1) = - x(6) * dx1theta - x(7) * dx2theta - x(8) * dx3theta - x(9) * dx4theta;
grad(2) = - x(8) * dx3alpha - x(9) * dx4alpha;
grad(3) = - x(8) * dx3dV - x(9) * dx4dV - x(10) * dx5dV;
grad(4:end) = gradU;

end