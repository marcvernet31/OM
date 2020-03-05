function [d,b] = uo_descent_direction(isd, icg, irc, nu, x_1, x, g)
%   Search Direction:
%       isd=1 : GM; isd=2 : CGM; isd=3 : BFGS; 
%       icg=1 : FR; icg=2 : PR+;
%       irc=0 : no restart; irc=1 : RC1; irc=2 : RC2;

    b=0;
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
            b= (g(x)' * (g(x) - g(x_1))) / norm(g(x_1))^2;
        end
        d = -g(x) + b*d_1;
   
        
    elseif(isd == 3)
        %QuasiNewton
    end
    
end

