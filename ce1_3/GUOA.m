%Method: gradient(=0), newton(=1)
function [opt,k] = GUOA(x,f,g, gg, almax,almin,rho,c1,c2,iW, max_iterations, method)
    k = 1;
    while norm(g(x)) >= 10^-6 && k <= max_iterations
        d = -g(x);
        if(method == 1)
            d = d / gg(x);
        end
    al = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);
    
    x_aux = x;
    x = x+al*d;k = k+1;
    fprintf('%d   %d\n', (x_aux+0.739085)/(x+0.739085), (x_aux+0.739085)/((x+0.739085)^2))
    end
    opt = x;
end


