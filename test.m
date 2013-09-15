% test

addpath('src');

close all;
clear all;
format long;

disp('Test:');
disp('1 - poprawnoœci modelu');
disp('2 - poprawnoœci wyliczania Hamiltonianu');
disp('3 - poprawnoœci wyliczania Psi');
disp('4 - poprawnoœci wyliczania gradientu');
n = input('Wybierz test: ');

if n == 1
    disp('testModelu');
    
    [param, ~] = parametry();
    
    thetaRad = 230 * 2 * pi / 360;
    dVE = 0.3;
    u = [0; 0; 0; 0; 0; 0];
    %u = [0.2; 0.2; 0.2; 0; 0; 0];
    
    zd = [thetaRad; dVE; u];
    T = 10;
    tau = [0 T/3 2*T/3 T];
    rho = 1;
    
    var = obliczenia(param.h0, tau);
    
    [ x, ~, ~, ~ ] = solver(zd, param, var, rho);
    wizualizacja(x);
elseif n == 2
    disp('testH');
    [param, ~] = parametry();
    
    x0 = [2 2 1 1 1 1 1 1 1];
    psi0 = [1 1 1 1 1 1 1 1 1];
    u = [1;1];
    
    ep=1e-6;
    
    [ dH, psi ] = testH(x0, psi0, u, param, ep);
    
    disp('Porównanie');
    disp ([dH, psi]);
    disp('Procentowo (%)');
    disp (abs((dH - psi) ./ psi) * 100);
elseif n == 3
    disp('testPsi');
    [param, ~] = parametry();
    
    x0 = [2 2 1 1 1 1 1 1 1];
    u = [0; 0; 0; 0];
    h0 = param.h0;
    T = 1;
    tau = [0 T/2 T];
    ep=1e-6;
    rho = 1;
    
    var = obliczenia(param.h0, tau);
    
    [ dQ, H ] = testPsi(x0, u, param, var, ep, rho);
    
    disp('Porównanie');
    disp ([dQ, H]);
    disp('Procentowo (%)');
    disp (abs((dQ - H) ./ H) * 100);
elseif n == 4
    disp('testGrad');
    [param, ~] = parametry();
    
    thetaRad = 230 * 2 * pi / 360;
    dVE = 0.3;
    u = [0; 0; 0; 0; 0; 0];
    %u = [0.2; 0.2; 0.2; 0; 0; 0];
    
    zd = [thetaRad; dVE; u];
    
    T = 10;
    tau = [0 T/3 2*T/3 T];
    ep=1e-6;
    rho = 1;
    
    var = obliczenia(param.h0, tau);
    
    [ dQ, grad ] = testGrad(zd, param, var, ep, rho);
    
    disp('Porównanie');
    disp ([dQ, grad]);
    disp('Procentowo (%)');
    disp (abs((dQ - grad) ./ grad) * 100);
else
    disp('B³¹d!');
end
