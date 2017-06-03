function rho_eff = density_eff(c,rho_i,rho_m)
%density_eff calculates eff density
%   c   = volumfraction of inclusion in matrix
%   rho = density
%   _i  = inclusion
%   _m  = matrix
rho_eff = c*rho_i + (1-c)*rho_m;
end