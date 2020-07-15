%Solució: alpha = 0.5; iWout = 2;
%CORRECTE

Q = [4 0; 0 1]; 
f  = @(x) (1/2)*x'*Q*x; 
g  = @(x) Q*x;x  = [1;2]; 
d = [-4;-2];
almax= 1.0; almin= 10^-6; 
rho = 0.5; c1 = 0.1; c2 = 0.5;
iW=1


[al,iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);


disp('alpha: ');
disp(al);
disp('iWout: ');
disp(iWout);
