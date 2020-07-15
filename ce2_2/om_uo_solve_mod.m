
function [xk,dk,alk,iWk,betak,Hk, check_p] = om_uo_solve_mod(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu)
    %FIRST DERIVATIVE METHOD
    %   Stopping criterium:
    %       epsG: conv. tolerance
    %       kmax: maxim. iterations
    %   BLS:
    %       iW=0 : ELS; iW=1 : WC; iW=2 : SWC;
    %   Search Direction:
    %       isd=1 : GM; isd=2 : CGM; isd=3 : BFGS; 
    %       icg=1 : FR; icg=2 : PR+;
    %       irc=0 : no restart; irc=1 : RC1; irc=2 : RC2; 

    %check_p: (només per GM)
%   1 si es compleix perpendicularitat entre d(k) i d(k+1)
%   0 si no

    I = [1,0; 0,1];
    H = I;
    xk = [x]; dk = []; alk = [];
    iWk = []; betak = []; Hk = [h]; check_p = [];
    k = 0; x_1 = x; d_1 = 0; p_k = 1;
    
    while norm(g(x)) >= epsG && k < kmax
        
        
        [d, b] = uo_descent_direction(isd, icg, irc, nu, x_1, x, g, d_1, H, k);
        [al, iWi] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);
        
        x_1 = x;
        x = x + al * d;
        
        %Actualitzacions per GM
        if(isd == 1 && k>1)
            if(dot(d, d_1) ~= 0)
                p_k = 0;
            end
        end
        %Actualitzacions pel Quasi-Newton
        if(isd==3)
            s = x - x_1;
            y = g(x) - g(x_1);
            p = 1 / (y'*s);
            H = (I - p*s*y') * H * (I - p*y*s') + p*s*s';
        end
        d_1 = d;
        xk = [xk, x]; dk = [dk, d]; alk = [alk, al];
        iWk = [iWk, iWi]; betak = [betak, b]; %Hk = [Hk, H];
        check_p = [check_p, p_k];
        k = k + 1;
    end
end
