%--------------------------------------------------------------------------
% Title: velocityprofile
% Author: Sergio Carrillo De Hert
% Affiliation: University College Dublin
% Last modified: March 05, 2024
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Description: Function that calculates the velocity for every rradial position 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Dependencies: 
%       N.A.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Usage: 
% Input data: 
%           Ri      Inner radius of the annulus
%           Ro      Outer radius of the annulus
%           v_mean  Mean axial velocity
%           r       Radial positions
% Output data:
%           velocity at radial position r
%--------------------------------------------------------------------------


function vz = velocityprofile(Ri, Ro, v_mean, r)
%Velocity profile for annular flow
    Kappa=Ri/Ro;                                    %Inner to outer radius ratio
    Lambda=((1-Kappa^2)/(2*log(Kappa^-1)))^0.5;     %Dimensionless radial position of maximum velocity
    vz=2*v_mean/(Kappa^2+1-2*Lambda^2)*(1-(r/Ro)^2+2*Lambda^2*log(r/Ro));
end