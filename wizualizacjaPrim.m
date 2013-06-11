function wizualizacjaPrim(x, t)
% rysowanie trajektorii po unieruchomieniu Ksiê¿yca

global mu

figure;
hold on;

cs = [cos(t); sin(t)]';
p = zeros(length(t),2);
p(:,1) = x(:,1) .* cs(:,1) + x(:,2) .* cs(:,2);
p(:,2) = - x(:,1) .* cs(:,2) + x(:,2) .* cs(:,1);

plot(p(:,1), p(:,2), '-r');
plot(p(1,1), p(1,2), 'xr');
plot(p(end,1), p(end,2), 'xr');
plot(-mu /(1 + mu),0, 'ob');
plot(1 / (1 + mu), 0, 'ob');
xlabel('x');
ylabel('y');
title('Trajektoria lotu w nieinercjalnym uk³adzie wspó³rzêdnych');

end