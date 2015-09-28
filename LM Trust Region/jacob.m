function [ Ja ] = jacob( fname,x)

h=1.e-6;
n=length(x);
Ja=zeros(length(feval(fname,x)),n);

for j = 1:length(x)
    e = zeros(n,1);
    e(j) = 1;
    fxh = feval(fname,x+h*e);
    fx = feval(fname,x);
    Ja(:,j) = (1/h)*(fxh-fx);
end
