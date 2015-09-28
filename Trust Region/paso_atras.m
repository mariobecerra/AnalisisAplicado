function [ alpha ] = paso_atras(fname,x, p, pend, pN)

% Arturo �vila 124990
% Emilio Morones 124311
% Guillermo Schulz 122375

% B�squeda de l�nea usando las condiciones de las rectas en el punto
% x con la direcci�n p.

% In
%   fname = cadena con el nombre de la funci�n en Matlab.
%   x = vector columna en R^n que representa el punto donde se desciende.
%   p = vector de descenso de opreden n.
%   pend = n�mero real con la derivada direccional de fname en x a lo largo de
%          del vector p.
% Out
%   alfa = n�mero real entre [0,1] donde x + alfa*p cumpla las
%          condiciones de las rectas.

c1 = 0.5;
c2 = 0.6;
jmax = 10;
jiter = 0;
alpha = 1;

    rx = feval(fname, x); 
    fx = (rx' * rx)/2; 
    
    xt = x + alpha*p;
    rxt = feval(fname, xt); 
    fxt = (rxt' * rxt)/2; 

while((fxt > fx + c1*alpha*pend || fxt < fx + c2*alpha*pend) && alpha >= 0.1 && jiter < jmax )
    alpha = (-pend*alpha^2)/(2*(fxt-fx-alpha*pend));
    xt = x + alpha*p;
    rxt = feval(fname, xt); 
    fxt = (rxt' * rxt)/2; 
    jiter = jiter + 1;
end

if jiter == jmax
    alpha = 0.1;
end













end

