function out = rhs(x, u)
% wyliczanie prawej strony rownania modelu dla RK4
global mu restmu C1

rE = norm([x(1) + mu, x(2)]); 
rM = norm([x(1) - restmu, x(2)]);

out = zeros(1,5);
out(1) = x(3);
out(2) = x(4);
out(3) = x(1) + 2 * x(4) - restmu * (x(1) + mu) / rE^3 - mu * (x(1) - restmu) / rM^3 + u(1) * cos(u(2)) / x(5);
out(4) = x(2) - 2 * x(3) - restmu * x(2) / rE^3 - mu * x(2) / rM^3 + u(1) * sin(u(2)) / x(5);
out(5) = - C1 * u(1);
end