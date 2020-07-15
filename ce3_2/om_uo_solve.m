
function [wo, niter] = om_uo_solve(x,f,g,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu,delta)
    %   Stopping criterium:
    %       epsG: conv. tolerance
    %       kmax: maxim. iterations
    %   BLS:
    %       iW=0 : ELS; iW=1 : WC; iW=2 : SWC;
    %   Search Direction:
    %       isd=1 : GM; isd=2 : CGM; isd=3 : BFGS(QN);  
    %
    %       icg=1 : FR; icg=2 : PR+;
    %       irc=0 : no restart; irc=1 : RC1; irc=2 : RC2; 
    %  !! H : matriu per quasi-Newton
    %  !! h : hessiana
    


    I = eye(35);
    H = I;
    xk = [x]; dk = []; alk = [];
    iWk = []; betak = []; al=0; iWi=0;
    k = 1; x_1 = x; d_1 = 0; tauk = [];
    Hk(:,:,1) = H;
    
    while norm(g(x)) >= epsG && k < kmax
        
        %Direcció de descens
        [d, b, Bk, tau] = uo_descent_direction(isd, icg, irc, nu, x_1, x, g, d_1, H, k, h,delta);
        
        %Longitud de pas
        [alpha,iout] = uo_BLSNW32(f,g,x,d,almax,c1,c2,kBLSmax,epsal)
        [al, iWi] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);
        
        x_1 = x;
        x = x + al * d;
        %Actualitzacions pel Quasi-Newton
        if(isd==3)
            s = x - x_1;
            y = g(x) - g(x_1);
            p = 1 / (y'*s);
            H = (I - p*s*y') * H * (I - p*y*s') + p*s*s';
            Hk(:,:,k) = H;
        end
        d_1 = d;
        iWk = [iWk, iWi];
        k = k + 1;
    end
    w = x;
    niter = k;
end
