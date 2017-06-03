function mu_eff = threePhaseModel_shear(c,nu_i,mu_i,nu_m,mu_m)
%threePhaseModel_shear solves A*(mu_eff/mu_m)^2 + 2B*(mu_eff/mu_m) + C = 0
%   c = volumefraction of inclusions in matrix. 
%   nu = poisson ratio
%   mu = shear modulus
%   _i = inclusion
%   _m = matrix

    n1 = ((49-50*nu_i*nu_m)*(mu_i/mu_m - 1)+35*(mu_i/mu_m)*(nu_i - 2*nu_m) ... 
         + 35*(2*nu_i-nu_m)); 
    n2 = 5*nu_i*(mu_i/mu_m - 8)+7*(mu_i/mu_m + 4);
    n3 = (mu_i/mu_m)*(8-10*nu_m)+(7-5*nu_m);
    
    A = (8*(mu_i/mu_m - 1)*(4-5*nu_m)*n1*c^(10.0/3)...
        - 2*(63*(mu_i/mu_m - 1)*n2+2*n1*n3)*c^(7.0/3)...
        + 252*(mu_i/mu_m - 1)*n2*c^(5.0/3)...
        - 50*(mu_i/mu_m - 1)*(7-12*nu_m+8*nu_m^2)*n2*c ...
        + 4*(7-10*nu_m)*n2*n3);
        
    B = (-2*(mu_i/mu_m - 1)*(1-5*nu_m)*n1*c^(10.0/3)...
        + 2*(63*(mu_i/mu_m - 1)*n2+2*n1*n3)*c^(7.0/3)...
        - 252*(mu_i/mu_m - 1)*n2*c^(5.00/3)...
        + 75*(mu_i/mu_m - 1)*(3-nu_m)*n2*nu_m*c... 
        + (3.0/2)*(15*nu_m-7)*n2*n3);

    C = (4*(mu_i/mu_m - 1)*(5*nu_m-7)*n1*c^(10.0/3)... 
        - 2*(63*(mu_i/mu_m - 1)*n2+2*n1*n3)*c^(7.0/3)... 
        + 252*(mu_i/mu_m - 1)*n2*c^(5.0/3)...
        + 25*(mu_i/mu_m - 1)*(nu_m^2-7)*n2*c... 
        - (7+5*nu_m)*n2*n3);     
    
    discremant = (2*B)^2 - 4*A*C;
    
    if discremant > 0
        %disp('Two solutions');
        X1 = (-(2*B)+sqrt(discremant))/(2*A);
        X2 = (-(2*B)-sqrt(discremant))/(2*A);
        
        X = max(X1,X2);
    
    elseif discremant == 0
        %disp('One solution');        
        X = (-(2*B)/(2*A));
    else
        %disp('This solution is imag');
        X = 0;
    end
    
    mu_eff = X*mu_m;
end

