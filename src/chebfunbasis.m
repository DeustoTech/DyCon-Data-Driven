function Y = chebfunbasis(X,n)
    
    if n == 0
        Y = 0*X + 1;
    elseif n==1
        Y = X;
    else
        Y = 2*X.*chebfunbasis(X,n-1) - chebfunbasis(X,n-2);
        
    end 
end