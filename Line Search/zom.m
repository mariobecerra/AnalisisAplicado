%%-------------------------------------------------------------------------
% FUNCI�N zom: Creamos la funci�n zoom como se expone en el libro Numerical
% Optimization de Nocedal y Wriht. Dado un intervalo (alfa_hi, alfa_lo), y 
% una funci�n encuentra un alfa_j intermedia tal que alfa_lo no viola la 
% condici�n de suficiente descenso(la que da el menor valor de la funci�n).
% 
% alfa_hi es elegida de tal modo que la derivada de la funci�n evaluada en
% alpha_lo por (alpha_hi - alpha_lo) sea menor que cero.
%                                                       "Salve esperanza"
%%-------------------------------------------------------------------------

function alpha_estrella = zom(x, p, alpha_lo, alpha_hi, fname, c1, c2)

%phi_prima = gradiente(fname, x);
g = gradiente(fname, x);
flag=1;

while flag
    alpha_lo
    alpha_j = 11*(alpha_lo + alpha_hi)/21;
    phi = feval(fname, x + p*alpha_j);
    
    if phi > feval(fname, x) + c1*alpha_j*g'*p || phi >= feval(fname, x + alpha_lo * p)
        alpha_hi = alpha_j
        
    else
        phip_alpha_j = gradiente(fname, x + alpha_j * p)'*p;
        
        if abs(phip_alpha_j) <= -c2* g'*p
            alpha_estrella=alpha_j;
            flag=0;
        end
        
        if phip_alpha_j * (alpha_hi - alpha_lo) >=0 
            alpha_hi = alpha_lo;
        end
        alpha_lo = alpha_j;
    end
end

end
    
            
          
      