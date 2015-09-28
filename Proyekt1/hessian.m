function [H] = hessian(fname,x)
% Aproximaci�n a la matriz hessiana por diferecencias
% hacia adelante.

% ITAM
% Zeferino Parada
% 21 de agosto de 2015
%In
% fname .- cadena con el nombre de la funcion.
% x .- vector columna en R^n donde se aproxima la hessiana.
% Out
% H .- matriz sim�trica de orden n.

n = length(x);         % dimensi�n del problema
fx = feval(fname,x);   % funci�n evluada en x
h = 1.e-06;            % tama�o de paso para las diferencias hacia adelante
H = zeros(n);          % espacio para la matriz hessiana

for i = 1:n
    for j = 1:i
        xt = x; xt(i) = xt(i) + h;
         H(i,j) = - feval(fname,xt);
         xt = x; xt(j) = xt(j) + h;
         H(i,j) = H(i,j) - feval(fname,xt);
         xt = x; xt(i) = xt(i) +h; xt(j) = xt(j) + h;
         H(i,j) = (H(i,j) + feval(fname,xt) + fx)/ (h^2);
         if ( j ~= i)  H(j,i) = H(i,j); end
     end
 end
 
 
 