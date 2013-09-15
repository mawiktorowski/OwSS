% g³ówny plik solvera

addpath('src');

close all;
clear all;
format long e;

thetaRad = 250;
thetaRad = degtorad(thetaRad);
dVE = 0.3;
n = 10; % iloœæ przedzia³ów
u = zeros(2*n,1);
zd0 = [thetaRad; dVE; u];
T = 10; % czas od którego przeszukiwanie zostaje rozpoczête
krokT = 10; % krok przeszukiwania
[param, opts] = parametry();
Qhist = [];
Thist = [];

r = 20; % d³ugoœæ kroku próbnego
alfa = 3;
Ne = 5; %iloœæ kroków ekspansji
beta = 0.5; % wspó³czynnik kontrakcji
epK = 2e-10;
Nz = 5; %iloœæ kroków z³otego podzia³u

% poszukiwanie rozwi¹zañ dopuszczalnych
% musi byæ przeszukiwane odpowiednio gêsto, aby nie przegapiæ minimum
while true
    tau = linspace(0,T,n+1);
    [ Q, zd, tau, dop ] = funkcjaKary(zd0, tau, param, opts);
    Thist = [Thist T];
    Qhist = [Qhist Q];
    if dop
        %disp('Przerwano');
        Td = T;
        Qd = Q;
        break;
    end
    T = T + krokT;
end

z1 = T;
Q1 = Q;
step = 1;
tp = (sqrt(5) - 1) / 2;

while true
    switch(step)
        case 1
            j = 1;
            z2 = z1 + r;
            tau = linspace(0,z2,n+1);
            [ Q2, zd, tau, dop ] = funkcjaKary(zd0, tau, param, opts);
            Thist = [Thist z2];
            Qhist = [Qhist Q2];
            if Q2 >= Q1
                step = 2;
            else
                step = 3;
            end
        case 2
            if z2 - z1 < epK
                % rozwi¹zanie znajduje siê w przedziale [z1, z2)
                fprintf('OSI¥GNIÊTO MAKSYMALN¥ DOK£ADNOŒÆ KONTRAKCJI\n');
                break;
            end
            Q3 = Q2;
            z2 = beta * z2;
            tau = linspace(0,z2,n+1);
            [ Q2, zd, tau, dop ] = funkcjaKary(zd0, tau, param, opts);
            Thist = [Thist z2];
            Qhist = [Qhist Q2];
            if Q2 >= Q1
                step = 2;
            else
                z3 = z2 / beta;
                step = 8;
            end
        case 3
            j = j + 1;
            if j > Ne
                Topt = z3;
                Qopt = Q3;
                tau = linspace(0,z3,n+1);
                [ Qopt, zdopt, tau, dop ] = funkcjaKary(zd0, tau, param, opts);
                step = 9;
            else
                z3 = alfa * z2;
                tau = linspace(0,z3,n+1);
                [ Q3, zd, tau, dop ] = funkcjaKary(zd0, tau, param, opts);
                Thist = [Thist z3];
                Qhist = [Qhist Q3];
                if Q3 >= Q2
                    j=0;
                    step = 5;
                else
                    step = 4;
                end
            end
        case 4
            Q1 = Q2;
            Q2 = Q3;
            z1 = z2;
            z2 = z3;
            step = 3;
        case 5
            j = j + 1;
            if j > Nz
                step = 8;
            else
                step = 6;
            end
        case 6
            if (z3 - z2) > (z2 - z1)
                step = 7;
            else
                zp = z3 + tp * (z2 - z3);
                tau = linspace(0,zp,n+1);
                [ Qp, zd, tau, dop ] = funkcjaKary(zd0, tau, param, opts);
                Thist = [Thist zp];
                Qhist = [Qhist Qp];
                if Qp >= Q2
                    Q3 = Qp;
                    z3 = zp;
                else
                    Q1 = Q2;
                    Q2 = Qp;
                    z1 = z2;
                    z2 = zp;
                end
                step = 5;
            end
        case 7
            zp = z1 + tp * (z2 - z1);
            tau = linspace(0,zp,n+1);
            [ Qp, zd, tau, dop ] = funkcjaKary(zd0, tau, param, opts);
            Thist = [Thist zp];
            Qhist = [Qhist Qp];
            if Qp >= Q2
                Q1 = Qp;
                z1 = zp;
            else
                Q3 = Q2;
                Q2 = Qp;
                z3 = z2;
                z2 = zp;
            end
            step = 5;
        case 8
            licz = (Q1 * (z2^2 - z3^2) + Q2 * (z3^2 - z1^2) + Q3 * (z1^2 - z2^2));
            mian = (Q1 * (z2 - z3) + Q2 * (z3 - z1) + Q3 * (z1 - z2));
            zm = 0.5 * licz / mian;
            tau = linspace(0,zm,n+1);
            [ Qm, zd, tau, dop ] = funkcjaKary(zd0, tau, param, opts);
            Thist = [Thist zm];
            Qhist = [Qhist Qm];
            if Qm < Q2
                Topt = zm;
                Qopt = Qm;
                zdopt = zd;
            else
                Topt = z2;
                tau = linspace(0,z2,n+1);
                [ Qopt, zdopt, tau, dop ] = funkcjaKary(zd0, tau, param, opts);
            end
            step = 9;
        case 9
            var = obliczenia(param.h0, tau);
            [ x, Q, ~, ~] = solver(zdopt, param, var, 0);
            wizualizacja(x);
            wyniki = [Qhist; Thist];
            [ ~, ind ] = sort(wyniki(2, :));
            wyniki = wyniki(:, ind);
            figure;
            hold on;
            plot(wyniki(2,:),wyniki(1,:), '-b');
            plot(wyniki(2,:),wyniki(1,:), 'xk');
            plot(Td, Qd, 'or');
            plot(Topt, Qopt, 'ok');
            hold off;
            xlabel('T(TU)');
            ylabel('Q');
            title('WskaŸnik Q (bez kary) w funkcji czasu');
            sterowanie(zdopt, tau, param);
            break;
    end
end
