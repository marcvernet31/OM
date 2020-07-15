%Solució: alpha = 0.0635; iWout = 3;
%CORRECTE

f = @(x) 40*sin(pi*x(1)/2+pi/2)*cos(pi*x(2)/4) + x(1)^2 + x(2)^2;
g = @(x) [ 20*pi*cos(pi*x(1)/2+pi/2)*cos(pi*x(2)/4) + 2*x(1);
    -10*pi*sin(pi*x(1)/2+pi/2)*sin(pi*x(2)/4) + 2*x(2)];
x  = [-3;0]; 
d = [15;0]; 
almax= 1.0; almin= 10^-6; 
rho = 0.5; c1 = 0.1; c2 = 0.5;
iW=2


[al,iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);


disp('alpha: ');
disp(al);
disp('iWout: ');
disp(iWout);