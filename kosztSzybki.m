function Q = kosztSzybki(x0, h0, tau, u)
% szybkie wyliczanie wskaznika jakosci bez zapamietywania trajektorii

dtau = diff(tau);
n = ceil(dtau./h0);
h = dtau./n;
h2 = h/2;
h3 = h/3;
h6 = h/6;
cn = cumsum([1 n]);

us = zeros(2,1);

x = x0;
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