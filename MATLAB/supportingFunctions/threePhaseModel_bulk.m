function k_eff = threePhaseModel_bulk(c,k_i,k_m,mu_m)
%threePhaseModel_bulk calculates the eff bulk modulus
%   c   = volumefraction of inclusions in matrix.
%   k   = bulk modulus 
%   mu  = shear modulus
%   _i  = inclusion
%   _m  = matrix
k_eff = k_m + ((k_i-k_m)*c/(1+(1-c)*((k_i-k_m)/(k_m+(4.0/3)*mu_m))));
end
