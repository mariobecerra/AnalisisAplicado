% ITAM
% An�lisis Aplicado con el doctor Zeferino Parada
% 21. august 2015
% Morales Mendoza, Rodrigo Andr�s 000124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ x, iter, W ] = MetodoNewton( fname, x )
% M�todo de m�ximo descenso con b�squeda de l�nea
n = length(x);
maxiter = 50;
maxjter = 20;
tol = 1.e-06;
iter = 0;
%jter = 0;
%c1 = 0.1;
c1 = 0.0001;
W = x; %Aqu� se ir�n guardando los puntos que genera la sucesi�n.
gx = gradiente( fname, x );
ng= norm(gx);
beta = 1e-3;% para la hessiana spd
maxith = 10;
while (ng > tol && iter < maxiter)
    H = hessian(fname,x);
    %%% Checar condicion de H.
% % % %     condition= cond(H);
% % % %     if condition>1e10
% % % %         disp('error en condicionamiento')
% % % %         break;
% % % %     end
    %--------------------------Fin check condicion H
    %Checar H spdf (posit def)
    aii=diag(H);
    a0=min(aii)
    if a0>0
        tk = 0;
        H
    else
        tk = -a0 + beta;
    end
    for k = 1:maxith
        [L p] = chol(H+tk*eye(n));
        if p<0
            H
            L'*L
            H = L'*L;
            break;
        else
            tk=max(2*tk,beta);
        end
    end
    if k == maxith
        disp('aqui se llego al maxith')
    end
    %---------------------Fin de checar si H spdf (posit def)
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
    %--------------------- fin de graficacion
    % B�squeda de l�nea
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Aqu'i ir'a lo de zoom etc pag 60 nocedal
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
    W = [W x];
    iter = iter + 1;
    gx = gradiente( fname, x );
    ng= norm(gx);
end

end