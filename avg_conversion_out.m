%--------------------------------------------------------------------------
% Title: Avg_conversion_out
% Author: Sergio Carrillo De Hert
% Affiliation: University College Dublin
% Last modified: March 05, 2024
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Description: Function that calculates the Q-weighted cross-sectional average conversion
% for the PDE model
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Dependencies: 
%           flowprofile.m
%           solvemasspde.m
%               masspde.m
%                   velocity profile.m
%                   LVPRAfunction.m
%               massic.m
%               massbc.m
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Usage: 
% Input data: Ri, Ro, v_mean, pDm, phik, N_LDF, kappa_PC, kappa_tot, L, rsteps, zsteps
%           Ri          Inner radius of the annulus
%           Ro          Outer radius of the annulus
%           v_mean      Mean axial velocity
%           pDm         -log base 10 of the molecular diffusivity
%           coefficient
%           phik        Product of quantum yield and kinetic constant
%           N_LDF       Flow rate of photons emitted by LDF
%           kappa_PC    Extinction coefficient of PC
%           kappa_tot   Total extinction coefficient
%           L           Length of reactor
%           rsteps      Discrete number of radial positions
%           zsteps      Discrite number of axial positions
% Output data:
%           TotX_final (Q) is the Q-weighted fractional conversion at the
%           outlet of the CAP-Flow system. Produces a vector of
%           size(v_mean)
%--------------------------------------------------------------------------

function TotX_final=avg_conversion_out(Ri, Ro, v_mean, pDm, phik, N_LDF, kappa_PC, kappa_tot, L, rsteps, zsteps)
    Dm=10^-(pDm);                                   %Obtain Dm
    q=flowprofile(Ri,Ro,v_mean,rsteps);             %Call flow rate function

    for j=1:length(v_mean)
        sol = solvemasspde(Ri, Ro, v_mean(j), Dm, phik, N_LDF, kappa_PC(j), kappa_tot(j), L, rsteps, zsteps); %Call solvemasspde
        F_P(j,:) = sol(end,:) .* q(j,:);                    %Conversion(r,z) multiplied by flow rate(r)
        TotX_final(j) = sum(F_P(j,:), 2) ./ sum(q(j,:),2);  %Divide former by total flow rate
    end
end