function wizualizacja(x)
% rysowanie trajektorii

global mu restmu

figure;
hold on;
plot(x(:,1),x(:,2), '-k');
plot(x(1,1),x(1,2), 'xk');
plot(x(end,1),x(end,2), 'xk');
plot(- mu,0, 'or');
plot(restmu, 0, 'or');
hold off;
xlabel('x(DU)');
ylabel('y(DU)');
title('Trajektoria lotu');

end