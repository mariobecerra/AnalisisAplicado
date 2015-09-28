function [x,iter] = Metodo_Descenso_Newton(fname,x)

% Encuentra minimos locales de funciones dos veces diferenciables
% de R^n a R donde el punto inicial esta muy cerca de un minimo local.
% La dirección de descenso es la dirección de Newton.

%In
%   fname = cadena con el nombre de la funcion.
%   x = vector columna de dimension n.

% Out
%   x = vector columna de dimension n, es la aproximacion al
%     minimo local.
%   iter = numero de iteraciones que tarda el metodo

% parámetros
tol = 1.e-08; % tolerancia para la norma del gradiente.
maxiter = 2000;    % número máximo de iteraciones externas permitidas
c1 = 1.e-04;
c2 = 0.9;

% valores iniciales
iter = 0;
g  = gradiente(fname,x);
normg    = norm(g);

disp(' iter       r           K           p0        f(x*)      ||f´(x*)|| ')
disp('----------------------------------------------------------------------')

while ( normg > tol && iter < maxiter)  % parte iterativa
    
    H = hessian(fname,x);  % H es simétrica
    m = min(eig(H));
    
    if m <= 0
        H = H +(abs(m)+1)*eye(length(H));    % H es s.p.d.
    end
    
    pN = H\-g;  % al ser s.p.d. Matlab utiliza Cholesky para resolver el sistema
    pN = pN/norm(pN);
    %pend = g'*pN;
    alpha  = line_search( fname, x, pN, c1, c2, 1);
    x = x + alpha*pN
    g = gradiente(fname,x);
    fx = feval(fname,x);
    normg = norm(g);
    iter = iter + 1;
    disp(sprintf('%3.0f \t %2.5f \t %2.5f \t %2.5f \t %2.5f \t %2.5f \t %2.5f',iter,x(1),x(2),x(3),fx,normg)  )
    
end

end