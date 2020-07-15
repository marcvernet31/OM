function [d,b] = uo_descent_direction(isd, icg, irc, nu, x_1, x, g, d_1, H, k)
%   Search Direction:
%       isd=1 : GM; isd=2 : CGM; isd=3 : BFGS; 
%       icg=1 : FR; icg=2 : PR+;
%       irc=0 : no restart; irc=1 : RC1; irc=2 : RC2;
    
    %Nombre de iteracions per RC1
    nRestart = 10;

    %Restart conditions
    RC1 = @(k) mod(k, nRestart) == 0; 
    RC2 = @(x) abs(g(x)'*g(x_1))/norm(g(x))^2 >= nu;

    b=0; I= [1,0; 0,1];
    
    %Restart
    if(irc == 1 && RC1(k))
        d = -g(x);
        %Falta contar les iteracions i canviar
    elseif(irc == 2 && RC2(x))
        d = -g(x);
        
    %No restart    
    else
        if(isd == 1)
            %GradientMethod
            d = -g(x);
        
        elseif(isd == 2)
            %ConjugatedGradientMethod
            if(icg == 1)
                %FR
                b = (g(x)' * g(x)) / norm(g(x_1))^2;
            elseif(icg == 2)
                %PR+
                b = max(0, (g(x)' * (g(x) - g(x_1))) / norm(g(x_1))^2);
            end
            d = -g(x) + b*d_1;
   
        elseif(isd == 3)
            %QuasiNewton
            d = -H*g(x);
        end
    end
end

