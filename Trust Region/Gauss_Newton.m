function [x, iter] = Gauss_Newton(fname,x0)
% Resuelve el problema de mínimos cuadrados con funciones no lineales
% utilizando el método de Gauss-Newton con la modificación de Levenberg-Marquadt. 
% Min f(x) = (1/2) * ( r(x)' * r(x) ) con r(x) función de Rn --> Rm
%
%In:
% fname.- nombre de función de residuales a minimizar
% x.- punto inicialt
% tol.- tolerancia del método
% maxiter.- número máximo de iteraciones.
%
%Out:
% x.- vector solución
% iter.- número total de iteraciones
%
% Consultado en 
% Nocedal, Jorge; Wright, Stephen J. (2006). Numerical Optimization, 2nd Edition. Springer. ISBN 0-387-30303-0,
% Capítulo 10: Least-Squares Problems
%
% ITAM, Septiembre 2015
% Análisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% José Carlos Castro 127049
% José Manuel Incera 125360
% Rodrigo Morales Mendoza 124341
%
% 

% parámetros iniciales
tol = 1.e-8;                   % tolerancia a la norma del gradiente
maxiter = 1000;                  % número máximo de iteraciones
iter = 0;                       % contador de iteraciones
% valores iniciales
x = x0;
rx = feval(fname,x);
Jx = jacob(fname,x);        % Matriz jacobiana
salida = norm(Jx'*rx);          % vector gradiente
n = length(x);

while( salida > tol && iter < maxiter)

    % Resuelva el problema de mínimos cuadrados
    %  Min (1/2) || J*p + rx ||_2^2

    % [p,k] = MiGradCon(Jx'*Jx,Jx'*(-rx)); %Minimizar la norma utilizando gradiente conjugado
    H = Jx'*Jx;

    b = eig(H); % eigenvalores de la aproximación de la Hessiana
    if (sum( b < 0 ) > 0) 
        % Si alguno de los eigenvalores es menor a cero, se le suma a B
        % en la diagonal el valor absoluto del mínimo eigenvalor para
        % que todos sean positivos. 
        % "two-dimensional subspace minimization", Nocedal cap. 4
        H = H + (1 + abs(min(b))) * eye(n);
    end

    
    b =  - Jx' * rx;
    p = H\b;

    t = 1.0;
    x = x + t*p;                   % nuevo punto
    rx = feval(fname,x);           % nuevo residual
    Jx = jacob(fname,x);       % nueva matriz jacobiana
    salida = norm(Jx'*rx);         % nuevo criterio de salida
    iter = iter +1;                % incremento en las iteraciones.
  
end



end