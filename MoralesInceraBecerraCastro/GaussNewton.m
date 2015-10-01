% ITAM, Septiembre 2015
% Análisis Aplicado con el doctor Zeferino Parada
% Mario Becerra 124362
% José Carlos Castro 127049
% José Manuel Incera 125360
% Rodrigo Andrés Morales Mendoza 124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proyecto: Curva de crecimiento logístico
% función de Newton
% iPad y Trigo
function [ x, iter ] = GaussNewton(fname, tol, maxiter, x)
% Resuelve el problema de mínimos cuadrados con funciones no lineales
% utilizando el método de Levenberg-Marquardt. 
% Min f(x) = (1/2) * ( r(x)' * r(x) ) con r(x) función de Rn --> Rm
%%%
%In:
% fname.- funktionName de residuales a minimizar
% x.- punto inicial
% tol.- tolerancia del método
% maxiter.- máximo de iteraciones.
%%%
%Out:
% x.- vector solución
% iter.- núm de iteraciones
% Consultado en 
% Nocedal, Jorge; Wright, Stephen J. (2006). Numerical Optimization, 2nd Edition. Springer. ISBN 0-387-30303-0,
% Capítulo 10: Least-Squares Problems
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
        %B = hessian('F',x);
        %eig(B) % Si se usa la Hessiana tiene un eigenvalor = -51000 y dos
        %pequeños
        %cond(B) % y además será mal condicionada
        B = Jx' * Jx; % Aprox de la Hessiana de f en x
        b = eig(B); % eigenvalores de la aprox. de la Hessiana
        if (sum( b < 0 ) > 0) 
            % Si alguno de los eigenvalores es menor a cero, se le suma a B
            % en la diagonal el valor absoluto del mín eigenvalor para
            % que todos sean positivos. 
            % "two-dimensional subspace minimization", Nocedal cap. 4
            B = B + (1 + abs(min(b))) * eye(n);
        end

        % Calcular dirccn de descenso
        while(flag == 0)
            p_temp = Doblez(B, grad, delta);

            r_x_p = feval(fname, x + p_temp);
            c = (rx' * rx) + c1 * grad' * p_temp;
            if( r_x_p' * r_x_p  <= c ) % Armijo (descenso suficiente)
               % Si se cumple la condcn de Armijo, toma p_temp como direccn
               p = p_temp;
               flag = 1;
            else
                % Si no se cumple se actualiza delta
                delta = max(delta/2, delta_lo);
                if(delta == delta_lo)
                    % Si la nueva delta coincide con delta_lo,
                    %       toma p_temp como dirccn
                    p = p_temp;
                    flag = 1;
                end
                    % Si la nueva delta no coincide con delta_lo, 
                    % calcula la
                    % funcn doblez de nuevo y se prueba la nieva dirccn
            end
        end
        flag = 0;
        

        % Actualizar radio de Trust Region (delta)
        mc = 0.5 * p' * B * p + grad' * p + fx; 
        % (modelo cuadr) para Funktn Trst Region
        f_x_p = feval(fname, x + p);
        rho = (f_x_p' * f_x_p - fx) / (mc - fx); 
        % Cociente de reduccn verdadera y reduccn predicha
        if (rho < 1.25 & rho > 0.75)
            delta = min(2 * delta, delta_hi);
        end

        % Actualizar variables
        x = x + p;    rx = feval(fname, x);    Jx = jacob(fname, x);
        grad = Jx' * rx;    fx = (rx' * rx) / 2;
        iter = iter+1;

    end
end

