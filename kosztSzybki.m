function Q = kosztSzybki(zd, h0, tau)
% szybkie wyliczanie wskaznika jakosci bez zapamietywania trajektorii

global C2 rE VE m mu

dtau = diff(tau);
n = ceil(dtau./h0);
h = dtau./n;
h2 = h/2;
h3 = h/3;
h6 = h/6;
cn = cumsum([1 n]);

x = zeros(1,5);
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
    f = @(x,u,t) rhs(x,u);
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

end