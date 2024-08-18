%--------------------------------------------------------------------------
% Title: LVPRAfunction
% Author: Sergio Carrillo De Hert
% Affiliation: University College Dublin
% Last modified: March 05, 2024
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Description: Function that calculates the LVPRA at r and z
% for the PDE model
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Dependencies: 
%           N.A.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Usage: 
% Input data: Ri, Ro, v_mean, pDm, phik, N_LDF, kappa_PC, kappa_tot, L, rsteps, zsteps
%           Ri          Inner radius of the annulus
%           kappa_tot   Total extinction coefficient
%           N_LDF       Flow rate of photons emitted by LDF
%           r           Radial position
%           z           Axial positions
% Output data:
%           LVPRA at r,z
%--------------------------------------------------------------------------

function LVPRA = LVPRAfunction( Ri, kappa_tot, N_LDF, r, z)
%LVPRA function
    LVPRA=(kappa_tot*N_LDF/(2*r*pi))*log(10)*10^(z-1)*exp(-kappa_tot*(r-Ri));
end