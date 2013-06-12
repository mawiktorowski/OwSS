function [ x, Q ] = BFGS(x0, h0, tau, u, ogr)
% algorytm BFGS

global ep0 ep1 ep2 epK0 epK1 MAX_ITER

x = u;

f = @(x) solverSzybki(x0, h0, tau, x);

epK = epK0;     % dok³adnoœæ kontrakcji - zale¿nie od kierunku d

R = 1;
iter = 1;

while(iter <= MAX_ITER)   
    [Q, grad] = f(x);
%     grad(grad == ogr(:,1) & grad > 0) = 0;
%     grad(grad == ogr(:,2) & grad < 0) = 0;
% mam wrazenie ze te znaki powinny byc odwrotnie
% stary kod prawdopodobnie do usuniecia
    %rzutowanie gradientu na ograniczenia
    for j = 1:length(x)
        if x(j) == ogr(j,1)
            grad(grad > 0) = 0;
        elseif x(j) == ogr(j,2)
            grad(grad < 0) = 0;
        end 
    end
    
    % ---- Krok 2 ----------------------------------
    if norm(grad) <= ep0 
        disp('KONIEC - MALA NORMA GRADIENTU. ');
        disp('ITERACJA: ');
        disp(iter);
        break 
    end
     
    % ---- Krok 3 -----------------------------------
    if R
        W = eye(length(x));
        epK = epK0;
    else
        r = grad - gradOld;
        s = x - xOld;
        Ws = W*s;
        W = W + (r*r')/(s'*r) - (Ws*s'*W)/(s'*Ws);
        epK = epK1;
    end
    
    d = W \ (-grad);
  
    % ---- Krok 4 -----------------------------------
    %if d'*grad > -min(ep1, ep2*norm(Q)^2)    
    if d'*grad > -min(ep1, ep2*norm(Q)^2)
        disp('ODNOWA ALGORYTMU');
        R = 1;      % idz do kroku 3
        iter = iter+1;
        continue    % przejœcie do kolejnej iteracji
    end
    
    % ---- Krok 5 -----------------------------------
    xOld = x;
    gradOld = grad;
    
    % ---- Krok 6 -----------------------------------
% line search !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%[x, Q] = lineSearch(x, d, Q, ogr, epK);
% ---- "rzutowanie wektora d" na ograniczenia
%d = max(d, ogr(:,1) - x);
%d = min(d, ogr(:,2) - x);
% stary kod prawdopodobnie do usuniecia
    for i = 1:length(x)
        d(x+d < ogr(i,1)) = ogr(i,1);
        d(x+d > ogr(i,2)) = ogr(i,2);
    end
lambda = 1;
% ---- wlasciwe poszukiwanie na kierunku
while(lambda > epK)
    xNew = x + lambda*d;
    Qnew = f(xNew);
    if(Qnew < Q)
        x = xNew;
        Q = Qnew;
        break
    else
        lambda = lambda/2;  % kontrakcja
    end
end 
% koniec line search !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    % ---- Krok 7 -----------------------------------
    if max(abs(xOld-x))   % je¿eli poprawa rozwi¹zania
        R = 0;     
    else 
       if  R; 
          disp('KONIEC - KIERUNEK NAJSZYBSZEGO SPADKU I BRAK POPRAWY');
          break; 
       else
           R = 1;   %idz do kroku 3
       end  
    end
    
    iter = iter+1
end
x
Q
end