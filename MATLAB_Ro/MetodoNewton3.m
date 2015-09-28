% ITAM
% Análisis Aplicado con el doctor Zeferino Parada
% 21. august 2015
% Morales Mendoza, Rodrigo Andrés 000124341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ x, iter, W ] = MetodoNewton( fname, x )
% Método de máximo descenso con búsqueda de línea
n = length(x);
maxiter = 50;
maxjter = 20;
tol = 1.e-06;
iter = 0;
%jter = 0;
%c1 = 0.1;
c1 = 0.0001;
W = x; %Aquí se irán guardando los puntos que genera la sucesión.
gx = gradiente( fname, x );
ng= norm(gx);
while (ng > tol && iter < maxiter)
    H = hessian(fname,x);
    %e=eig(H)   %un eig es muy negativo, tipo -50mil
    %%% Checar condicion de H.
% % % %     condition= cond(H);
% % % %     if condition>1e10
% % % %         disp('error en condicionamiento')
% % % %         break;
% % % %     end
    %--------------------------Fin check condicion H
    %Checar H spdf (posit def)
    H = hazspd(H);
    %e=eig(H); %checar que H es spdf
    %---------------------Fin de checar si H spdf (posit def)
    %%%%%%%%%%%%%%%
    %p = -gx; %%%%%%Ésta es la única condición que cambió
    p=-H\gx; %Dirección de Newton
    %%%%%%%%%%%%%%%
    %Graficación de fname a lo largo de x+p
    t = linspace(0,0.01,50)';
    ft=zeros(50,1);
    for k=1:50
        ft(k) = feval(fname, x+t(k)*p);
    end
    figure(1)
    plot(t,ft,'--b','LineWidth',3)
    title('Gráfica de búsqueda de línea')
    pause(0.5);
    %--------------------- fin de graficacion
    % Búsqueda de línea
    
    
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
    %---------------------------fin búsqueda de línea
    x = x + alfa*p;
    W = [W x];
    iter = iter + 1;
    gx = gradiente( fname, x );
    ng= norm(gx);
end

end