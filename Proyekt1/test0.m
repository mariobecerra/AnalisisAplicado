% for k=1:5
%     if(k==3)
%         break
%     end
% end
% if k==5
%     disp('si fue 5')
% end

%Este script nos permite ver si una matriz es SPD, sim pos def
clear
clc
close all
a =[5 2 -3; 2 5 4; -3 4 3];

[m n] = size(a);
[V, D]=eig(a); % a*V = V*D, con a simetrica V^-1 = V'  %tons a = V*D*V'
delta = 1e-2;
d=diag(D);
ti = zeros(n,1);
ti(d<delta)=delta-d(d<delta);
anew= V*(D+diag(ti))*V';

% Checking
[V,D]=eig(anew);
diag(D)

% % % %ee = min(e);
% % % %diag(a);
% % % %a = a+(abs(ee)+1)*eye(n);    % H es s.p.d.
% % % %e = eig(a)
% % % [L p] = chol(a+0*eye(n));
