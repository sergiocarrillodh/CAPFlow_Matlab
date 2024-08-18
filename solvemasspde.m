%--------------------------------------------------------------------------
% Title: solvemasspde
% Author: Sergio Carrillo De Hert
% Affiliation: University College Dublin
% Last modified: March 05, 2024
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Description: Solves PDE equation according to handle function masspde
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Dependencies: 
%       masspde.m
%           velocity profile.m
%           LVPRAfunction.m
%       massic.m
%       massbc.m
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Usage: 
% Input data: 
%           Ri          Inner radius of the annulus
%           Ro          Outer radius of the annulus
%           v_mean      Mean axial velocity
%           Dm          Molecular diffusivity coefficient
%           phik        Product of quantum yield and kinetic constant
%           N_LDF       Flow rate of photons emitted by LDF
%           kappa_PC    Extinction coefficient of PC
%           kappa_tot   Total extinction coefficient
%           L           Length of the reactor
%           rsteps      Discrete number of radial positions
%           zsteps      Discrete number of axial positions
% Output data:
%           Fractional conversion for each r,z position
%--------------------------------------------------------------------------



function solconvprofile = solvemasspde(Ri, Ro, v_mean, Dm, phik, N_LDF, kappa_PC, kappa_tot, L, rsteps, zsteps)
    r = linspace(Ri, Ro, rsteps);
    z = linspace(0, L, zsteps);
    solconvprofile = pdepe(1, @(r, z, u, dudr) masspde(r, z, u, dudr, Ri, Ro, v_mean, Dm, phik, N_LDF, kappa_PC, kappa_tot), @massic, @massbc, r, z); %Call pdepe function
end