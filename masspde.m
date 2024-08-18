%--------------------------------------------------------------------------
% Title: masspde
% Author: Sergio Carrillo De Hert
% Affiliation: University College Dublin
% Last modified: March 05, 2024
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Description: function handle that defines the coefficients of the PDE for
% each radial and axial position r,z
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Dependencies: 
%           velocity profile.m
%           LVPRAfunction.m
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Usage: 
% Input data: Ri, Ro, v_mean, phik, N_LDF, kappa_PC, kappa_tot, L, rsteps, zsteps
%           r           radial position
%           z           axial position
%           Ri          Inner radius of the annulus
%           Ro          Outer radius of the annulus
%           v_mean      Mean axial velocity
%           Dm          Molecular diffusivity coefficient
%           phik        Product of quantum yield and kinetic constant
%           N_LDF       Flow rate of photons emitted by LDF
%           kappa_PC    Extinction coefficient of PC
%           kappa_tot   Total extinction coefficient
% Output data:
%           Function handle vector
%--------------------------------------------------------------------------

function [c,f,s] = masspde(r,z,u,dudr, Ri, Ro, v_mean, Dm, phik, N_LDF, kappa_PC, kappa_tot)
    c = velocityprofile(Ri, Ro, v_mean, r);             %Call Velocity profile function 
    f=Dm.*dudr;                                         %Flux term in PDE
    %Call LVPRA function
    LVPRA= LVPRAfunction( Ri, kappa_tot, N_LDF, r, z);  %Call LVRPAfunction
    s=LVPRA.*phik.*(1-u);                               %Source term for first order kinetics
end