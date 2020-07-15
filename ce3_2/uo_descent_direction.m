function [d,b, Bk, tau] = uo_descent_direction(isd, icg, irc, nu, x_1, x, g, d_1, H, k, h, delta)
    %   Search Direction:
    %       isd=1 : GM; isd=2 : CGM; isd=3 : BFGS; isd=4 : NM; 
    %       isd=5 : MNM-SD; isd=6 : MNM-CMI
    %
    %       icg=1 : FR; icg=2 : PR+;
    %       irc=0 : no restart; irc=1 : RC1; irc=2 : RC2; 
    

    %Restart conditions
    nRestart = 10; %Nombre de iteracions per RC1
    RC1 = @(k) mod(k, nRestart) == 0; 
    RC2 = @(x) abs(g(x)'*g(x_1))/norm(g(x))^2 >= nu;

    b=0; Bk=[]; I=[1,0; 0,1]; tau = 0;
    
    %GM
    if(isd==1)
        d = -g(x);
        
    %CGM
    elseif(isd==2)
        %Restart
        if((irc == 1 && RC1(k)) || (irc == 2 && RC2(x)))
            d = -g(x);
        %No restart
        else
            if(icg == 1)%FR
                b = (g(x)' * g(x)) / norm(g(x_1))^2;
            elseif(icg == 2)%PR+
                b = max(0, (g(x)' * (g(x) - g(x_1))) / norm(g(x_1))^2);
            end
            d = -g(x) + b*d_1;
        end
        
    %BFGS(Q-N)
    elseif(isd==3)
        d = -H*g(x);
        
    %NM
    elseif(isd==4)
        d = -inv(h(x)) * g(x);
        
    %MNM-SD
    elseif(isd==5)
        [V, D] = eig(h(x));
        for i = 1:length(h(x));
            D(i,i) = max (D(i,i), delta);
        end
        Bk = V * D * V';
        d = -inv(Bk) * g(x);
  
    %MNM-CMI
    elseif(isd==6)
        laUB = norm(h(x), 'fro');
        it = 0; p=1; tau=0;
        while p ~= 0
            tau = (1.01 - 1/(2^it)) * laUB;
            [R,p] = chol(h(x) + tau*eye(2));
            it = it + 1;
        end
        Bk = R'*R;
        d = -inv(Bk) * g(x);
    end
end

