function [ Q, grad ] = solverSzybki(x0, h0, tau, u)
% solver szybki modelu rozwiazujacy rownania bez zapamietywania trajektorii
% zwraca wskaznik jakosci i gradient

global restmu K2 rM VM mr

dtau = diff(tau);
n = ceil(dtau./h0);
h = dtau./n;
h2 = h/2;
h3 = h/3;
h6 = h/6;
cn = cumsum([1 n]);

x = x0;
t = 0;

grad = zeros(2*length(dtau),1);
us = zeros(2,1);

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
    grad(j) = x(11);
    grad(length(dtau)+j) = x(12);
    x(11:12) = [0 0]; % gamma(t_i+1) = 0
end
end