% ITAM
% An�lisis Aplicado con el doctor Zeferino Parada
% 21. august 2015
% Morales Mendoza, Rodrigo Andr�s 000124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ x, iter, W ] = MetodoNewton( fname, x )
% M�todo de m�ximo descenso con b�squeda de l�nea

maxiter = 50;
maxjter = 20;
tol = 1.e-06;
iter = 0;
%jter = 0;
%c1 = 0.1;
c1 = 0.0001;
W = [x]; %Aqu� se ir�n guardando los puntos que genera la sucesi�n.
gx = gradiente( fname, x );
ng= norm(gx);
H = hessian(fname,x);
while (ng > tol && iter < maxiter)
    %%%%%%%%%%%%%%%
    %p = -gx; %%%%%%�sta es la �nica condici�n que cambi�
     
    p=-H\gx; %Direcci�n de Newton
    %%%%%%%%%%%%%%%
    %Graficaci�n de fname a lo largo de x+p
    t = linspace(0,0.01,50)';
    ft=zeros(50,1);
    for k=1:50
        ft(k) = feval(fname, x+t(k)*p);
    end
    figure(1)
    plot(t,ft,'--b','LineWidth',3)
    title('Gr�fica de b�squeda de l�nea')
    pause(0.5);
    % B�squeda de l�nea
    alfa = 1.0;
    xt = x+p;
    fx = feval(fname,x);
    gx = gradiente( fname, x );
    gx(2)
    pend = gx'*p;
    fxt = feval(fname,xt);
    jter = 0;
    while( fxt > fx+alfa*c1*pend && jter < maxjter)
        alfa = alfa*3/4;
        xt = x+alfa*p;
        fxt = feval(fname,xt);
        jter = jter+1;
    end
    %---------------------------fin b�squeda de l�nea
    x = x + alfa*p;
    H = hessian(fname,x);
    condition= cond(H)
    if condition>1e10
        break;
    end
    W = [W x];
    iter = iter + 1;
    gx = gradiente( fname, x );
    ng= norm(gx);
end

end