% Problem
clear;
a = -1;
f  = @(x) x(1)^2 + x(2)^4 + (a/2)*x(2)^2; 
g  = @(x) [ 2*x(1); 4*x(2)^3 + a*x(2)];
h  = @(x) [ 2 , 0; 0 , 12*x(2)^2 + a];
x1 = [2; 0.1*sqrt(abs(a)/12)];
%x1 = [0.5;1];

    %   Stopping criterium:
    %       epsG: conv. tolerance
    %       kmax: maxim. iterations
    %   BLS:
    %       iW=0 : ELS; iW=1 : WC; iW=2 : SWC;
    %   Search Direction:
    %       isd=1 : GM; isd=2 : CGM; isd=3 : BFGS; isd=4 : NM; 
    %       isd=5 : MNM-SD; isd=6 : MNM-CMI
    %
    %       icg=1 : FR; icg=2 : PR+;
    %       irc=0 : no restart; irc=1 : RC1; irc=2 : RC2; 
    %  !! H : matriu per quasi-Newton
    %  !! h : hessiana


% Stopping criterium:
epsG = sqrt(eps); kmax= 100;
% Linesearch:
almax= 1.0; almin= 10^-3; rho=0.5; c1=0.01; c2=0.9; iW= 1; 
% Searchdirection:
isd = 5; icg= 1; irc= 2 ; nu= 0.1; delta = 0.1;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimization
[xk,dk,alk,iWk,betak,Hk,tauk] = om_uo_solve (x1,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu,delta);
save('uo_SDM_CE31.mat','f','g','h'  ,'epsG','kmax','almax','almin','rho','c1','c2','iW','isd','icg','irc','nu','delta','xk','dk','alk','iWk','betak','Hk','tauk');

% Output
niter = size(xk,2); xo = xk(:,niter);
fk  = []; gk = []; rk = []; gdk = []; rk = []; Mk = []; la1k = zeros(1,niter); kappak = zeros(1,niter);
for k = 1:niter
    x = xk(:,k); fk = [fk,f(x)]; gk = [gk,g(x)]; if isd > 3 la1k(k) = min(eig(h(x))); end
    if isd == 4 & all(eig(h(x))>0) kappak(k) = cond(h(xk(:,k))); end
    if k < niter
        rk = [ rk,norm(g(xk(:,k+1))) / norm(g(  xk(:,k))) ]; Mk =[  Mk, norm(g(xk(:,k+1))) / norm(g(xk(:,k)))^2 ];
        gdk = [gdk,gk(:,k)'*dk(:,k)];
        if isd == 3 la1k=min(eig(Hk(:,:,k))); end
        if  isd >= 3 & isd ~=4  kappak(k) = cond(Hk(:,:,k)); end
    end
end
if isd <= 4  tauk(1:niter-1) = 0; end
if isd == 5  tauk(1:niter-1) = delta; end
if isd == 3  la1k(niter) = 0; end

fprintf('[uo_SDM_CE31]\n');
fprintf('   f= %s\n', func2str(f));
fprintf('   epsG= %3.1e, kmax= %4d\n', epsG,kmax);
if isd ~= 4 fprintf('   almax= %2d, almin= %3.1e, rho= %4.2f, c1= %3.2f, c2= %3.2f, iW= %1d\n',almax,almin,rho,c1,c2,iW); 
end
fprintf('   isd= %1d\n',isd);
if isd == 2 fprintf('   icg= %1d, irc= %1d, nu= %3.1f\n',icg,irc,nu); end 
if isd == 5 fprintf('   delta= %3.1d\n',delta); end 
fprintf('   x1 = [ %+3.1e , %+3.1e ]\n', x1(1), x1(2));
fprintf('      k     g''*d       al iWla(1) del./tau    kappa    ||g||        f        r        M\n');
for k = 1:niter-1
    fprintf(' %6d %+3.1e %+3.1e  %1d %+3.1e %+3.1e %+3.1e %+3.1e %+3.1e %+3.1e %+3.1e\n', k, gdk(k), alk(k), iWk(k), la1k(k), tauk(k), kappak(k), norm(gk(:,k)), fk(k), rk(k),Mk(k));
end
fprintf(' %6d                      %+3.1e                   %+3.1e %+3.1e\n', niter, la1k(niter), norm(gk(:,niter)), fk(niter));
fprintf('      k     g''*d       al iWla(1) del./tau    kappa    ||g||        f        r        M\n');
fprintf('   x* = [ %+3.1e , %+3.1e ]\n', xo(1), xo(2));
fprintf('[uo_SDM_CE31]\n');
xylim = [0 0 0 0];
subplot(2,2,1); uo_solve_plot(f, xk, gk, xylim, 1, 18); subplot(2,2,2); uo_solve_plot(f, xk, gk, xylim, 2, 18);
subplot(2,2,3); plot(rk(1:niter-1),'-o'  ); title('r^k'); subplot(2,2,4); plot(Mk(1:niter-1),'-x'); title('M^k');
