function [ x, iter ] = reg_conf(fname, tol, maxiter, x)


    % VVecindad de confianza
    delta = 10;
    deltaMax = 50; 
    eta = 0.1; 


    n = length(x);
    Jx = jacob(fname, x); %Jacobiana de r evaluada en el punto inicial
    rx = feval(fname, x); % Vector de residuales evaluado en el punto inicial
    grad = Jx' * rx; % Gradiente de f en el punto inicial
    fx = (rx' * rx)/2; % Valor de f en el punto inicial
    B = Jx' * Jx;

    iter = 0;
    
    while (norm(grad) > tol  && iter < maxiter) % ciclo del algoritmo
        
        m = min(eig(B)); 
        if m < 1.e-04 
            B = B +( abs(m) + 0.01 )*eye(n); % B es s.p.d.
        end
        
        p = Doblez( B, grad, delta ); 
        
        r_x_p = feval(fname, x + p);
        
        mc = 0.5 * p' * B * p + grad' * p + fx;
        rho = (r_x_p' * r_x_p - fx) / (mc - fx);

        if rho < 0.25    % actualiza radio de búsqueda
            delta = 1/4 * delta;
        else 
            if (rho > 0.75  &&  norm(p) == delta)
                delta = min( 2*delta, deltaMax );
            end
        end

        if rho > eta % actualizas todo
            disp('actualiza')
            p
            x = x + p;
            rx = feval(fname, x);    
            Jx = jacob(fname, x);    
            grad = Jx' * rx;    
            fx = (rx' * rx) / 2; 
            B = Jx' * Jx;
        end

        iter = iter + 1;

    end

end















