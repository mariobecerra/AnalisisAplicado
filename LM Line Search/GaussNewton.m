function [x, iter] = GaussNewton(fname ,tol, maxiter,x)
% Resuelve el problema de mínimos cuadrados con funciones no lineales
% utilizando el método de Gauss-Newton con la modificación de Levenberg-Marquadt
% usando búsqueda de línea
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


    rx = feval(fname,x); % Vector de residuales evaluado en x
    Jx = jacob(fname,x); %Jacobiana de r evaluada en x
    norma_grad = norm(Jx'*rx); % norma del vector gradiente de f evaluado en x
    n = length(x);

    iter = 0;

    while( norma_grad > tol && iter < maxiter)

        % Resuelva el problema de mínimos cuadrados
        %  Min (1/2) || J*p + rx ||_2^2

        % [p,k] = MiGradCon(Jx'*Jx,Jx'*(-rx)); %Minimizar la norma utilizando gradiente conjugado
        H = Jx'*Jx;

        eigen = eig(H); % eigenvalores de la aproximación de la Hessiana
        if (sum( eigen < 0 ) > 0) 
            % Si alguno de los eigenvalores es menor a cero, se le suma a B
            % en la diagonal el valor absoluto del mínimo eigenvalor para
            % que todos sean positivos. 
            % "two-dimensional subspace minimization", Nocedal cap. 4
            H = H + (1 + abs(min(eigen))) * eye(n);
        end


        grad = Jx' * rx;
        p = - H\grad;
        
        % Búsqueda de línea
        alpha = 2;
        c = 1e-4;
        jmax = 10;
        jiter = 0;
        fx = (rx' * rx)/2;
        xt = x + alpha*p;
        rxt = feval(fname, xt);
        fxt = (rxt' * rxt)/2;
        pend = grad'*p;
        while( (fxt > fx + c*alpha*pend | alpha >= 1 ) && jiter < jmax )
            alpha = alpha/2;
            xt = x + alpha*p;
            rxt = feval(fname, xt);
            fxt = (rxt' * rxt)/2;
            jiter = jiter + 1;
        end

        %Nuevo punto
        x = x + alpha*p;

        %Evaluar parámetros en nuevo punto
        rx = feval(fname,x); % Residuales
        Jx = jacob(fname,x); % Jacobiana de r
        norma_grad = norm(Jx'*rx); % norma del gradiente de f
        iter = iter +1;

    end



end