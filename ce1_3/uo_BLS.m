function [al,iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW)
    WC1 = @(al) f(x+al*d) <= f(x) + c1*al*g(x)'*d;
    WC2 = @(al)g(x+al*d)'*d >= c2*g(x)'*d;
    SWC2 = @(al) abs(g(x+al*d)'*d) <= c2*abs(g(x)'*d);
    WC = @(al) WC1(al) & WC2(al);
    SWC = @(al) WC1(al) & SWC2(al);
    
    al = almax;
    if(iW == 1)
        
        while(not(WC(al)) && al > almin)
            al = al*rho;
        end
        
    elseif(iW == 2)
        
        while(not(SWC(al)) && al > almin)
            al = al*rho;
        end
        
    end

 iWout = SWC(al) + WC(al) + WC1(al);
end
