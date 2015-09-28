function [ p ] = Doblez2( pN, pC,Delta )
% Función de doblez para región de confianza usando
% el punto de Cauchy y la dirección de Newton

 z1 = norm(pN);
 z2 = norm(pC);
 
 if (z1 <= Delta)
      p = pN;
     % disp('Newton');
 else
     if (z2 >= Delta)
         p = Delta*(pC/z2);
         %disp('Cauchy');
     else
         A = (pN-pC)'*(pN-pC);
         B = 2*pC'*(pN-pC)  ;
         C = z2^2 - Delta^2  ;
         r = roots([A B C]);
         if (r(1) >0)
             ts = r(1);
         else
             ts = r(2);
         end
         p = pC + ts*(pN-pC);
         %disp('Doblez');
     end
 end
 
end

