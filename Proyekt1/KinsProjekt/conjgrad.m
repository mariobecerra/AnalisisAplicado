function [x] = conjgrad(A,b,x)

    r=b-A*x;
    p=r;
    r1=r'*r;
 
    for i=1:10000000
        Ap=A*p;
        a=r1/(p'*Ap);
        x=x+a*p;
        r=r-a*Ap;
        r2=r'*r;
        if sqrt(r2)<1e-10
              break;
        end
        p=r+r2/r1*p;
        r1=r2;
    end

end

