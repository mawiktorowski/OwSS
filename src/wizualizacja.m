function wizualizacja(x)
% rysowanie trajektorii

figure;
hold on;
%axis equal;
plot(x(:,1),x(:,2), '-.r');
plot(x(1,1),x(1,2), 'xr');
plot(x(end,1),x(end,2), 'xr');
plot(x(:,3),x(:,4), '-k');
plot(x(1,3),x(1,4), 'xk');
plot(x(end,3),x(end,4), 'xk');
plot(0,0, 'or');
hold off;
xlabel('x(DU)');
ylabel('y(DU)');
title('Trajektoria lotu');

end
