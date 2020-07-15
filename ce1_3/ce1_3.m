%Solució: x* = -0.739085
%CORRECTE

f = @(x) x^2 + 2 * sin(x);
g = @(x) 2 * cos(x) + 2 * x;
gg = @(x) 2 - 2 * sin(x);
almax = 1; almin = 0.01; rho = 0.75;
c1 = 0.1; c2 = 0.5;

max_iterations = 100;
method = 1;
iW = 2;
%Punt inicial pot ser qualsevol?
x = 5;

[opt, k] = GUOA(x,f,g, gg, almax, almin, rho, c1, c2, iW, max_iterations, method);

disp('iteracions ');
disp(k);
disp('optim: ')
disp(opt);