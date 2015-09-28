function [ alpha ] = line_search( fname, x, p, c1, c2, alpha_max )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

alpha_lo = 0;
alpha_hi = rand() * alpha_max;
i = 1;
g = gradiente(fname, x);
flag = 1;

while flag
    phi = feval(fname, x + alpha_hi * p);
    if (phi > feval(fname, x) + c1 * alpha_hi * g'*p || ( phi >= feval(fname, x + alpha_lo * p) && i>1))
        flag = 0;
        alpha = zom(x, p, alpha_lo, alpha_hi, fname, c1, c2);
    end
    phi_prima = gradiente(fname, x + alpha_hi * p)'*p;
    if ( abs(phi_prima) <= -c2 * g'*p )
        flag = 0;
        alpha = alpha_hi;
    end
    if (phi_prima >= 0)
        flag = 0;
        alpha = zom(x, p, alpha_hi, alpha_lo, fname, c1, c2);
    end
    temp = alpha_lo;
    alpha_lo = alpha_hi;
    alpha_hi = (alpha_max - temp) * rand() + temp; %Escoge un número entre temp y alpha_max
    i = i+1;
end


end

