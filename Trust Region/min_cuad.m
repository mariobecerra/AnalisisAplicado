function [ x, iter ] = min_cuad(fname, tol, maxiter, x)
% Resuelve el problema de mínimos cuadrados con funciones no lineales
% utilizando el método de Levenberg-Marquardt. 
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

    delta_lo = 1.e-05; delta_hi = 1; c1=1.e-04;


    delta = delta_hi;
    n = length(x);
    Jx = jacob(fname, x); %Jacobiana de r evaluada en el punto inicial
    rx = feval(fname, x); % Vector de residuales evaluado en el punto inicial
    grad = Jx' * rx; % Gradiente de f en el punto inicial
    fx = (rx' * rx)/2; % Valor de f en el punto inicial

    iter = 0;
    flag = 0;

    while(norm(grad) > tol && iter < maxiter)

        B = Jx' * Jx; % Aproximación de la Hessiana de f en x
        b = eig(B); % eigenvalores de la aproximación de la Hessiana
        if (sum( b < 0 ) > 0) 
            % Si alguno de los eigenvalores es menor a cero, se le suma a B
            % en la diagonal el valor absoluto del mínimo eigenvalor para
            % que todos sean positivos. 
            % "two-dimensional subspace minimization", Nocedal cap. 4
            B = B + (1 + abs(min(b))) * eye(n);
        end

        % Calcular dirección de descenso
        while(flag == 0)
            p_temp = Doblez(B, grad, delta);

            r_x_p_temp = feval(fname, x + p_temp);
            rx_temp = feval(fname,x);
            c = (rx_temp' * rx_temp) + c1 * grad' * p_temp;
            if( r_x_p_temp' * r_x_p_temp  <= c ) %Condición de Armijo (descenso suficiente)
               % Si se cumple la condición de Armijo, toma p_temp como dirección
               p = p_temp;
               flag = 1;
            else
                % Si no se cumple la condición de Armijo, se actualiza delta
                delta = max(delta/2, delta_lo);
                if(delta == delta_lo)
                    % Si la nueva delta coincide con delta_lo, toma p_temp como dirección
                    p = p_temp;
                    flag = 1;
                end
                    % Si la nueva delta no coincide con delta_lo, calcula la
                    % función doblez de nuevo y se prueba la nieva dirección
            end
        end
        flag = 0;
        

        % Actualizar radio de región de confianza (delta)
        mc = 0.5 * p' * B * p + grad' * p + fx; % Función modelo de región de confianza (modelo cuadrático)
        f_x_p = feval(fname, x + p);
        rho = (f_x_p' * f_x_p - fx) / (mc - fx); % Cociente de reducción verdadera y reducción predicha
        if (rho < 1.25 & rho > 0.75)
            delta = min(2 * delta, delta_hi);
        end

        % Actualizar variables
        x = x + p;    rx = feval(fname, x);    Jx = jacob(fname, x);    grad = Jx' * rx;    fx = (rx' * rx) / 2;
        iter = iter+1;

    end
end

