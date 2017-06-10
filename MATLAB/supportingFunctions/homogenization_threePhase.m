function [composite] = homogenization_threePhase(core,shell, matrix, c_sphere, c_composite) 
% homogenization_threePhase used the Three Phase Method (TPM) to create
% a homogeneous representation of a heterogeneous composite. Returns a
% material struct representing the homogeneous represeted composite.
%
% core:         Material struct with core.rho, core.nu and core.E
%               rho: Density, nu: Poisson's ratio, E: Young's modulus.
% shell:        Material struct.
% matrix:       Material struct.
% c_sphere:     Volume fraction of the core in the particle.
% c_composte:   Volume fraction of the particle in the composite.

%% Functions
shear   = @(E,nu)(E/(2*(1+nu)));
bulk    = @(E,nu)(E/(3*(1-2*nu)));
poisson = @(k,mu)((3*k-2*mu)/(2*(3*k+mu)));
youngs  = @(k,mu)(9*k*mu/(3*k+mu));
longVel = @(k,mu,rho)(sqrt((k+(4.0/3)*mu)/rho));

%% Material calculations
shell.mu        = shear(shell.E,shell.nu);
shell.k         = bulk(shell.E,shell.nu);
core.k          = bulk(core.E,core.nu);
core.mu         = shear(core.E,core.nu);
matrix.mu        = shear(matrix.E,matrix.nu);
matrix.k         = bulk(matrix.E,matrix.nu);

%% Homogenization of particle
sphere.mu       = threePhaseModel_shear(c_sphere,core.nu,core.mu,shell.nu,shell.mu);
sphere.k        = threePhaseModel_bulk(c_sphere,core.k,shell.k,shell.mu);
sphere.rho      = density_eff(c_sphere,core.rho,shell.rho);
sphere.nu       = poisson(sphere.k,sphere.mu);
sphere.cL       = longVel(sphere.k,sphere.mu,sphere.rho);
sphere.imp		= sphere.cL*sphere.rho;


%% Homogenization of composite
composite.mu    = threePhaseModel_shear(c_composite,sphere.nu,sphere.mu,matrix.nu,matrix.mu);
composite.k     = threePhaseModel_bulk(c_composite,sphere.k,matrix.k,matrix.mu);
composite.rho   = density_eff(c_composite,sphere.rho,matrix.rho);
composite.cL    = longVel(composite.k,composite.mu,composite.rho);
composite.imp   = composite.cL*composite.rho;
