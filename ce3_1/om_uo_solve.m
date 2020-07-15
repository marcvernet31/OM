
function [xk,dk,alk,iWk,betak,Hk, tauk] = om_uo_solve(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu,delta)
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
    


    I = [1,0; 0,1];
    H = I;
    xk = [x]; dk = []; alk = [];
    iWk = []; betak = []; al=0; iWi=0;
    k = 1; x_1 = x; d_1 = 0; tauk = [];
    Hk(:,:,1) = H;
    
    while norm(g(x)) >= epsG && k < kmax
        
        %Direcció de descens
        [d, b, Bk, tau] = uo_descent_direction(isd, icg, irc, nu, x_1, x, g, d_1, H, k, h,delta);
        
        %Longitud de pas
        if(isd==4)
            al = 1;
        else
            [al, iWi] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);
        end
        
        x_1 = x;
        x = x + al * d;
        %Actualitzacions pel Quasi-Newton
        if(isd==3)
            s = x - x_1;
            y = g(x) - g(x_1);
            p = 1 / (y'*s);
            H = (I - p*s*y') * H * (I - p*y*s') + p*s*s';
            Hk(:,:,k) = H;
        %Actualitzacions per MNM-SD i MNM-CMI
        elseif(isd==5 || isd==6)
            Hk(:,:,k) = Bk;
        end
        d_1 = d;
        xk = [xk, x]; dk = [dk, d]; alk = [alk, al];
        iWk = [iWk, iWi]; betak = [betak, b]; tauk = [tauk, tau];
        k = k + 1;
    end
end
