function wizualizacja(t, x, param)
% rysowanie trajektorii

figure;
hold on;
axis equal;
cs = [param.D * cos(param.omega * t) param.D * sin(param.omega * t)];
plot(cs(:,1),cs(:,2), '-.r');
plot(cs(1,1),cs(1,2), 'xr');
plot(cs(end,1),cs(end,2), 'xr');
plot(x(:,1),x(:,2), '-k');
plot(x(1,1),x(1,2), 'xk');
plot(x(end,1),x(end,2), 'xk');
plot(0,0, 'or');
hold off;
xlabel('x(DU)');
ylabel('y(DU)');
title('Trajektoria lotu');

end
