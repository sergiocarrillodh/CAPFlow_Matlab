%--------------------------------------------------------------------------
% Title: CAP-Flow system PDE Model: Effect of Q and Dm on X_A
% Author: Sergio Carrillo De Hert
% Affiliation: University College Dublin
% Last modified: March 05, 2024
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Description: The aim of this script analyses the effect of varying the
% volumetric flow rate ($Q$) and the molecular diffusion coefficient on
% conversion ($D_\text{m}$). The plug flow model is also simulated  for
% comparison.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Dependencies: 
%       avg_conversion_out.m
%       avg_conversion_outzero
%       solvemasspde.m
%       flowprofile.m
%       masspde.m
%       massbc.m
%       massic.m
%       velocityprofile.m
%       LVPRAfunction.m
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Usage: 
% The input data include: 
%           1. C_PC to be tested
%           2. Q range to be tested
%           3. Dm range to be tested
%           3. Geometry of the CAP-Flow system
%           4. Photon absorption properties of matrix and PC
%           5. Photon flow rate (by actinometry)
%           6. Kinetic contants (from
%           Fitter)
% The output includes:
%           1. Plot of X(Dm) for varying Q for PDE and plug-flow models
%--------------------------------------------------------------------------

clc;
clear;

%--------------------------------------------------------------------------
%Concentration of A  and photocatalyst
%--------------------------------------------------------------------------
CA0=0.4;                                            %Concentration of A [mol/L]
CP0=CA0*(0.07/100);                                 %Initial concentration of Ru(bpy)3(PF6)2 [mol/L]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Flow rates to be simulated in [mL/min]
%--------------------------------------------------------------------------
Q_mlmin=[0.5, 0.75, 1, 2, 4];
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Tubing dimensions
%--------------------------------------------------------------------------
Di_in=1/8;                                          %FEP outer diameter in [in]
Do_in=3/8-0.035*2;                                  %SS internal diameter in [in]
L=0.94;                                             %Length of LDF between inlet and outlet of CAP-Flow [m]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Photon Absorption coefficients
%--------------------------------------------------------------------------
alpha_PC=3.32*10^6;                                 %Naperian Molar absorptivity of photocatalyst [L/(mol m)]
kappa_matrix=110.5;                                 %Naperian extintion coefficient of the matrix without photocatalyst [m^-1]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Actinometry results, photon flow rate.
%--------------------------------------------------------------------------
N_Acti=7.7324e-7;                                   %Photon flow rate by Actinometry [einstein/s]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Molecular diffusion coefficient and kinetic constants (from Fitter.m)
%--------------------------------------------------------------------------
load('constants.mat', 'phik');
pDm=linspace(7,12,20);                              %Molecular diffusion coefficient in [m^2/s]
% phik=0.0168;                                        %quantum yield times kinetic constant
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%End of input data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Dimensions-conversions and derived geometric variables
%--------------------------------------------------------------------------
Ri_in=Di_in/2;                                      %FEP outer radius in [in]
Ri=Ri_in*0.0254;                                    %FEP outer radius in [m]
Ro_in=Do_in/2;                                      %SS internal radius in [in]
Ro=Ro_in*0.0254;                                    %SS internal radius in [m]
A_cross=pi*(Ro^2-Ri^2);                             %Cross sectional area in [m^2]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Naperian extintion coefficient
%--------------------------------------------------------------------------
kappa_PC=alpha_PC.*CP0;                             %Naperian Extintion coefficient of PC [m^-1]
kappa_tot=kappa_PC+kappa_matrix;                    %Napierian Extintion coefficient of mixture [m^-1]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Actinometry results, photon flow rate.
%--------------------------------------------------------------------------
N_LDF=N_Acti/(10^-0.07-(10^-1));                    %Photon flow rate by whole LDF [einstein/s]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Mesh for PDE
%--------------------------------------------------------------------------
rsteps=100;                                         %number of steps in the radial direction
zsteps=50;                                          %number of steps in the axial direction
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Flow rate conversion to SI units, mean velocity and mean residence time
%--------------------------------------------------------------------------
Q=(Q_mlmin.*10^-6)./60;                             %Flow rate in [m^3/s]
v_mean=Q/A_cross;                                   %mean velocity in [m/s]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%PDE model
%--------------------------------------------------------------------------
for j=1:length(Q)
   for i=1:length(pDm)
        TotX_final(j, i)=avg_conversion_out(Ri, Ro, v_mean(j), pDm(i), phik, N_LDF, kappa_PC, kappa_tot, L, rsteps, zsteps);
        Dm(j,i)=10^-pDm(i);
   end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Plug-flow CAP-Flow model
%--------------------------------------------------------------------------
for j=1:length(Q)
    for i=1:length(pDm)
        TotX_final_mean(j,i)=1-(exp((phik*N_LDF/(Q(j)))*(1)*(-1+exp(-kappa_tot*(Ro-Ri)))*(10^-(1-L)-10^-1)));
    end
end

%--------------------------------------------------------------------------
%Purely convective model
%--------------------------------------------------------------------------
for j=1:length(Q)
        TotX_final_zero(j)=avg_conversion_outzero(Ri, Ro, v_mean(j), phik, N_LDF, kappa_PC, kappa_tot, L, rsteps, zsteps);
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Figures
%--------------------------------------------------------------------------

repeated_Tot_X_zero = repmat(TotX_final_zero', 1, size(Dm, 2));

% Define a color map for the curves
colors = lines(length(Q));

figure(1);
hold on;
% Get handle to current axes and turn on the box around the plot
ax = gca; ax.Box = 'on';
for i = 1:length(Q)
    color = colors(i, :);
    h1=semilogx(Dm(i,:), TotX_final(i, :), '-', 'Color', color);
    h2=semilogx(Dm(i,:), TotX_final_mean(i, :), '--', 'Color', color);
    %h3=semilogx(Dm(i,:), repeated_Tot_X_zero(i,:), ':', 'Color', color);
    % Add legend label only for the continuous line
    legendLabels{i} = sprintf('%0.2f', Q_mlmin(i));
    % Exclude the dashed line from the legend
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    %set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
hold off;
xlabel('$D_m$ [m$^2$ s$^{-1}$]', 'Interpreter', 'latex');
ylabel('$\left\langle X_A \right\rangle$ [-]', 'Interpreter', 'latex');
lgd=legend(legendLabels, 'Interpreter', 'latex');
% Add a title to the legend
title(lgd, '$Q$ [mL min$^{-1}$]', 'Interpreter', 'latex');
%legend('0.50 mL min$^{-1}$', '0.75 mL min$^{-1}$', '1.00 mL min$^{-1}$', '2.00 mL min$^{-1}$', '4.00 mL min$^{-1}$', '6.00 mL min$^{-1}$', 'Interpreter', 'latex' );