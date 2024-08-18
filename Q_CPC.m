%--------------------------------------------------------------------------
% Title: CAP-Flow system PDE Model: Effect of Q and C_PC on X_A
% Author: Sergio Carrillo De Hert
% Affiliation: University College Dublin
% Last modified: March 05, 2024
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Description: The aim of this script is to analyze the effect of varying
% the feed volumetric flow rate and the photocatalyst concentration on
% conversion and on the coefficient of variation (CV) at the outlet of the
% CAP-Flow system
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Dependencies: 
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
%           2. Q to be tested
%           3. Geometry of the CAP-Flow system
%           4. Photon absorption properties of matrix and PC
%           5. Photon flow rate (by actinometry)
%           6. Molecular diffusion coefficient and kinetic contants (from
%           Fitter)
% The output includes:
%           1. Plot of X(C_PC) for varying Q
%           2. Coefficient of Variation (CV) as function of C_PC for
%           varying Q
%--------------------------------------------------------------------------

clc;
clear;

%--------------------------------------------------------------------------
%Concentration of A and simulated C_PC in percentage photon equivalents
%--------------------------------------------------------------------------
CA0=0.4;                                                                                %Concentration of A [mol/L]
CP0_eqperc=[0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1];  %Eq. percentage of PC
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Flow rates to be simulated in [mL/min]
%--------------------------------------------------------------------------
Q_mlmin=[0.5, 0.75, 1, 2, 4];
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Tubing dimensions
%--------------------------------------------------------------------------
Di_in=1/8;                                      %FEP outer diameter in [in]
Do_in=3/8-0.035*2;                              %SS internal diameter in [in]
L=0.94;                                         %Length of LDF between inlet and outlet of CAP-Flow [m]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Photon Absorption coefficients
%--------------------------------------------------------------------------
alpha_PC=3.32*10^6;                             %Naperian Molar absorptivity of photocatalyst [L/(mol m)]
kappa_matrix=110.5;                             %Naperian extintion coefficient of the matrix without photocatalyst [m^-1]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Actinometry results, photon flow rate.
%--------------------------------------------------------------------------
N_Acti=7.7324e-7;                               %Photon flow rate by Actinometry [einstein/s]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Molecular diffusion coefficient and kinetic constants (from Fitter.m)
%--------------------------------------------------------------------------
load('constants.mat', 'pDm', 'phik');
Dm=10^-pDm;                                 %Molecular diffusion coefficient in [m^2/s]
% phik=0.0168;                                    %quantum yield times kinetic constant
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
V_r=pi*(Ro^2-Ri^2)*L;                               %Volume reactor in [m^3]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Molar concentration of C_PC and naperian extintion coefficient
%--------------------------------------------------------------------------
CP0=CP0_eqperc.*(CA0/100);                          %Initial concentration of PC [mol/L]
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
r=linspace(Ri,Ro,rsteps);                           %Create linear space from Ri to Ro 
z=[linspace(0,L,zsteps)];                           %Create linear space from 0 to L
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Flow rate conversion to SI units, mean velocity and mean residence time
%--------------------------------------------------------------------------
Q=(Q_mlmin.*10^-6)./60;                             %Flow rate in [m^3/s]
v_mean=Q/A_cross;                                   %mean velocity in [m/s]
t_res=V_r./Q;                                       %Dimensionless time [-]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Calculate volumetric flow rate profile fraction for each r-step
%--------------------------------------------------------------------------
q=flowprofile(Ri,Ro,v_mean,rsteps);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Solve PDE model for each Q and C_PC
%--------------------------------------------------------------------------
for j=1:length(Q)
   for i=1:length(CP0_eqperc)
        sol = pdepe(1, @(r, z, u, dudr) masspde(r, z, u, dudr, Ri, Ro, v_mean(j), Dm, phik, N_LDF, kappa_PC(i), kappa_tot(i)), @massic, @massbc, r, z);
        X_outlet_r(:, j, i) = sol(end,:);
   end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Calculate CV for the CAP-Flow system's outlet for the different Q and C_PC
%--------------------------------------------------------------------------
for j=1:length(Q)
   for i=1:length(CP0_eqperc)
       TotX_final(j, i) = sum(X_outlet_r(:, j, i).*q(j,:)', 1) / Q(j);
       CoV_final(j,i)=((sum(q(j,:)'.*(X_outlet_r(:, j, i)-TotX_final(j, i)).^2,1)./(((rsteps-1)/rsteps)*Q(j)))^0.5)./TotX_final(j, i);
   end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Figures
%--------------------------------------------------------------------------
figure(1);
hold on; 
% Loop through Q values and plot TotX_final for each C_PC
for i = 1:length(t_res)
    plot(CP0_eqperc, TotX_final(i, :));
end
hold off;  % Disable plot overlay
xlabel('$C_{PC}$ [$\%$]', 'Interpreter', 'latex');
ylabel('$$\left\langle X_A \right\rangle$$ [-]', 'Interpreter', 'latex');
legend('0.50 mL min$^{-1}$', '0.75 mL min$^{-1}$', '1.00 mL min$^{-1}$', '2.00 mL min$^{-1}$', '4.00 mL min$^{-1}$', '6.00 mL min$^{-1}$', 'Interpreter', 'latex' ); 

figure(2);
hold on;
% Get handle to current axes
ax = gca;

% Turn on the box around the plot
ax.Box = 'on';
% Loop through Q values and plot TotX_final for each C_PC
for i = 1:length(Q_mlmin)
    plot(CP0_eqperc/100,CoV_final(i,:));
end
    xlabel('$C_{PC}$ [equiv.]', 'Interpreter', 'latex');
    ylabel('$CV$ [-]', 'Interpreter', 'latex');
    lgd=legend('0.50', '0.75', '1.00', '2.00', '4.00', 'Interpreter', 'latex' );
    % Add a title to the legend
title(lgd, '$Q$ [mL min$^{-1}$]', 'Interpreter', 'latex');