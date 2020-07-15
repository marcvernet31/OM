%Solució: alpha = 0.25; iWout = 2;
%CORRECTE

f = @(x) x(1)*exp(-x(1)^2-x(2)^2);
g = @(x) exp(-x(1)^2-x(2)^2)*[(1 -2*x(1)^2); (-2*x(1)*x(2))];
x  = [-1;1]; 
d = [3;-3]; 
almax= 1.0; almin= 10^-6; 
rho = 0.5; c1 = 0.1; c2 = 0.5;
iW=1


[al,iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);


disp('alpha: ');
disp(al);
disp('iWout: ');
disp(iWout);