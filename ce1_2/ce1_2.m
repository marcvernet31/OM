clear;
Q = [4 0; 0 1];
f = @(x) (1/2)*x'*Q*x; g = @(x) Q*x;
x = [1;2]; k = 1; iW =2;
almax = 1; almin = 0.01; rho = 0.5;
c1 = 0.1; c2 = 0.5;
b = [0 0; 0 0];

k = 1;
while norm(g(x)) >= 10^-6 && k <= 100
    d = -g(x);
    
    %Exact line search:
    %al = -((Q*x)' * d)/(d'*Q*d);
    
    %BLS with Wolfe Conditions:
    %al = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,1);
    
    %BLS with Strong Wolfe Conditions:
    al = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,2);
    
    fprintf('%3d %10.6f %10.6f %10.6f \n', k, x(1), x(2),al);
    x=x+al*d; k=k+1;
end
fprintf('%3d %10.6f %10.6f\n', k, x(1), x(2));


% El cas de Exact linial search fan falta 24 iteracions per arribar a [0,
% 0], i el resultat no és del tot precís.
% Amb Wolfe Conditions fan falta 4 iteracions, amb resultat precís.
%Amb Strong Wolfe Conditions fan falta 3 iteracions, amb resultat precís.
% En aquest cas, utilitzar exact line search no ens proporciona cap
% avantatge o millora.