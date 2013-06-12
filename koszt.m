function [Q, kara] = koszt(xT, T)
% wyliczanie wskaznik jakosci

global restmu K1 K2 rM VM mr

   function y1 = k1
       y1 = 1/4 * ((xT(1) - restmu)^2 + xT(2)^2 - rM^2)^2;    
   end 

   function y2 = k2
       y2 = 1/4 * ((xT(3) - xT(2))^2 + (xT(4) + xT(1))^2 - VM^2)^2;
   end 

   function y3 = k3
       y3 = 1/2 * ((xT(1) - restmu) * (xT(3) - xT(2)) + xT(2) * (xT(4) + xT(1)))^2;
   end

    function y4 = k4
        if xT(5) < mr
            y4 = 1/2 * (xT(5) - mr)^2;
        else
            y4 = 0;
        end
    end

kara = k1 + k2 + k3 + k4;
Q = - xT(5) + K1 * T + K2 * kara;
%[xT(5) T k1 k2 k3 k4]
end
