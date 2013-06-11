function out = rhsBack(x, u)
% wyliczanie prawej strony rownan w postaci kanonicznej
global mu restmu C1

rE = norm([x(1) + mu, x(2)]); 
rM = norm([x(1) - restmu, x(2)]); 

p31 = -1 + restmu * (rE^2 - 3 * (x(1) + mu)^2) / rE^5 + mu * (rM^2 - 3 * (x(1) - restmu)^2) / rM^5;
p41 = - 3 * restmu * (x(1) + mu) * x(2) / rE^5 - 3 * mu * (x(1) - restmu) * x(2) / rM^5;
p32 = p41;
p42 = -1 + restmu * (rE^2 - 3 * x(2)^2) / rE^5 + mu * (rM^2 - 3 * x(2)^2) / rM^5;

out = zeros(1,10);
out(1) = x(3);
out(2) = x(4);
out(3) = x(1) + 2 * x(4) - restmu * (x(1) + mu) / rE^3 - mu * (x(1) - restmu) / rM^3 + u(1) * cos(u(2)) / x(5);
out(4) = x(2) - 2 * x(3) - restmu * x(2) / rE^3 - mu * x(2) / rM^3 + u(1) * sin(u(2)) / x(5);
out(5) = - C1 * u(1);
% rownania sprzezone
out(6) = x(8) * p31 + x(9) * p41;
out(7) = x(8) * p32 + x(9) * p42;
out(8) = - x(6) + 2 * x(9);
out(9) = - x(7) - 2 * x(8);
out(10) = x(8) * u(1) * cos(u(2)) / x(5)^2 + x(9) * u(1) * sin(u(2)) / x(5)^2;
% dH/du
out(11) = out(8) * cos(u(2)) / x(5)^2 + out(9) * sin(u(2)) / x(5)^2;
out(12) = - out(8) * u(1) * sin(u(2)) / x(5)^2 + out(9) * u(1) * cos(u(2)) / x(5)^2;
end