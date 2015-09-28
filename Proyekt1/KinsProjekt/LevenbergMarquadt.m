function [ x,iter,fx,norma ] = LevenbergMarquadt(fname2,x)
% M�todo de minimizaci�n de funciones por algoritmo de regi�n de confianza.
% Recibe nombre de funci�n y punto inicial.

delta_max=1;
delta_min=1.e-05;
tol=1.e-06;
maxiter=100;
iter=0;
bandera=0;
delta=delta_max;
c1=1.e-04;
n = length(x);
Jx=jacobiana2(fname2,x);
rx=feval(fname2,x);
grad=Jx'*rx;
fx=.5*(rx'*rx);



while(norm(grad)>tol&&iter<maxiter)
    
    B=Jx'*Jx;
    b = eig(B);
    if (sum(b<0) > 0)
        B = B + (1+abs(min(b)))*eye(n);
    end
    
%     L = chol(B,'lower');
%     b1 = trin(L,-grad);
%     pN = tris(L',b1);
    pN= B\-grad;
    pC= -grad*(grad'*grad)/(grad'*B*grad);
    
    %Aceptaci�n de paso
    while(bandera==0)
        p2=Doblez2(pN,pC,delta);
        
        
        if(((feval(fname2,x+p2)'*feval(fname2,x+p2))<=(feval(fname2,x)'*feval(fname2,x))+c1*grad'*p2))
           p=p2;
           bandera=1;
        else
            delta=max(delta_min,delta*0.5);
            if(delta==delta_min)
                p=p2;
                bandera=1;
            end
        end
    end
    
    % Actualizaci�n del radio
    mcp=0.5*p'*B*p+grad'*p+fx;
    ro=((feval(fname2,x+p)'*feval(fname2,x+p))-fx)/(mcp-fx);
    
    if (ro > 0.75 && ro < 1.25)
        delta=min(2*delta,delta_max);
    end
    
    % Actualizaci�n de variables
    x=x+p;
    iter =iter+1;
    bandera=0;
    
    Jx=jacobiana2(fname2,x);
    rx=feval(fname2,x);
    grad=Jx'*rx;
    fx=.5*(rx'*rx);
    
end
norma=norm(grad);
end

