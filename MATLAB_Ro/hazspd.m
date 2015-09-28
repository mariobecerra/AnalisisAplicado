function [ anew ] = hazspd( a )
%Recibe matriz y la hace spd
[m n] = size(a);
[V, D]=eig(a); % a*V = V*D, con a simetrica V^-1 = V'  %tons a = V*D*V'
d=diag(D);
delta = 1e-2;
if min(d)>delta
    anew = a;
else
    ti = zeros(n,1);
    ti(d<0)=delta-d(d<0)
    %ti(d<delta)=delta-d(d<delta);
    anew= V*(D+diag(ti))*V';
end

end

