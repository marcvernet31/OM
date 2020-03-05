clear;
% Problem
f = @(x) x(1)^2 + x(2)^3 + 3*x(1)*x(2);
g = @(x) [ 2*x(1)+3*x(2); 3*x(2)^2 + 3*x(1)];
h = @(x) []; x = [-3;-1];
% Input parameters.
 % Stopping criterium:
epsG = 10^-6; kmax = 1500;
 % Linesearch:
almax = 2; almin = 10^-3; rho=0.5; c1=0.01; c2=0.45; iW = 2;
 % Search direction:
isd = 1; icg = 2; irc = 0 ; nu = 0.1;

% Optimization
[xk,dk,alk,iWk,betak,Hk] = om_uo_solve(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu);
% save('uo_FDM_CE21.mat','f','g','h','epsG','kmax','almax','almin','rho','c1','c2','iW','isd','icg','irc','nu','xk','dk','a lk','iWk','betak');

% Output
niter = size(xk,2); xo = xk(:,niter); fk = []; gk = []; rk = []; gdk = [];
for k = 1:niter
    fk = [fk,f(xk(:,k))]; gk=[gk,g(xk(:,k))]; % f(xk) and g(xk)
end
for k = 1:niter-1
    rk = [rk,(fk(k+1)-f(xo))/(fk(k)-f(xo))]; % Rate of convergence
end
for k = 1:niter-1
    gdk = [gdk,gk(:,k)'*dk(:,k)]; % Descent condition
end
rk=[rk,NaN]; gdk = [gdk,0];
fprintf('[om_uo_FDM_CE21]\n');
fprintf('epsG= %3.1e, kmax= %4d\n', epsG,kmax);
fprintf('almax= %2d, almin= %3.1e, rho= %4.2f\n',almax,almin,rho);
fprintf('c1= %3.2f, c2= %3.2f, iW= %1d\n',c1,c2,iW);
fprintf('isd= %1d, icg= %1d, irc= %1d, nu= %3.1f\n',isd,icg,irc,nu);
fprintf('k x(1) x(2) iW g''*d ||g|| r\n');
for k = 1:niter-1
    fprintf('%5d %7.4f %7.4f %3d %+3.1e %4.2e %3.1e\n', k, xk(1,k), xk(2,k), iWk(k), gdk(k), norm(gk(:,k)), rk(k));
end
fprintf(' k x(1) x(2) iW g''*d ||g|| r\n[om_uo_FDM_CE21]\n');
xylim=[0 0 0 0];
subplot(2,1,1); uo_solve_plot(f, xk, gk, xylim, 1, 0); subplot(2,1,2); uo_solve_plot(f, xk, gk, xylim, 2, 0);