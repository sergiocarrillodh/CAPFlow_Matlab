%--------------------------------------------------------------------------
% Title: flowprofile
% Author: Sergio Carrillo De Hert
% Affiliation: University College Dublin
% Last modified: March 05, 2024
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Description: Function that calculates the Q(i) for every r_i to r_{i-1} position for
% the PDE model
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
%           rsteps  Discrete number of radial positions
% Output data:
%           Flow rate (Q) with dimensions jxi, where j is given by the
%           number of v_mean and, i by the number of rsteps
%--------------------------------------------------------------------------

function q=flowprofile(Ri,Ro,v_mean,rsteps)
    %Discretized flow rate
    r = linspace(Ri, Ro, rsteps);                   %create vector with r(i) positions
    Kappa=Ri/Ro;                                    %inner-to-outer radius ratio
    Lambda=((1-Kappa^2)/(2*log(Kappa^-1)))^0.5;     %Dimensionless radial position of maximum velocity
    q=zeros(length(v_mean),rsteps);                 %Preallocate variables for speed

    for i = 2:rsteps                                % Start from 2, as we skip the first index 
        term1 = r(i)^2 - r(i - 1)^2 - (r(i)^4 - r(i - 1)^4) / (2*Ro^2 );
        term2 = r(i)^2 * (2*log(r(i)/Ro) - 1) - r(i - 1)^2 * (2*log(r(i - 1)/Ro) - 1);
        term3 = (2*pi / (Kappa^2 + 1 - 2 * Lambda^2)) * (term1 + Lambda^2 * term2);

        for j = 1:length(v_mean)
            q(j, i) = v_mean(j) * term3;
        end
    
    end

end