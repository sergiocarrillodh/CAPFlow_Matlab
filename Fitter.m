%--------------------------------------------------------------------------
% Title: CAP-Flow system PDE Model: Fitter
% Author: Sergio Carrillo De Hert
% Affiliation: University College Dublin
% Last modified: July 12, 2024
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Description: The aim of this script is to analyse the experimental data.
% It calculates the External ($Q$Y). The absorbed photon equivalents
% ($\eta_\text{a,eq}$). Obtains the molecular diffusivity coefficient 
% ($D_\text{m}$) and kinetic constants ($\phi_\lambda\cdot k$) from the
% experimental data. The parity plot and the experimental and modelled 
% $X_\text{A}$ are compared. Calculates $R^2$ for each data set.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Dependencies: 
%       avg_conversion_out.m
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
%       1. Experimental results, the photocatalyst concentration, flow rate and fractional conversion as input
%       2. Dimensions of CAP-Flow system for different geoemtries to the
%       ones speciefied in the manuscript
%       3. Photon flow rate obtained from actinometry
% The output includes:
%       A. Figures
%           1. Plot of X(Q)
%           2. Parity plot
%           3. Plot of X(Q); experimental vs predicted
%           4. Plot of X(absorbed photon equivalent)
%       B. Numerical - stored as variables
%           1. pDm and phik: -log10(molecular diffusivity) and kinetic
%           constants from experimental data using the LSE
%           2. Estimates external quantum yield
%           2. Calculate R^2 for each data set
%--------------------------------------------------------------------------

clc;
clear;

%--------------------------------------------------------------------------
%Experimental data
%First column: photocatalyst concentration in percentage equivalents
% Second column: flow rate in mL/min
% Third column: fractional conversion
%--------------------------------------------------------------------------
data0025=[
%    0.0025, 6.00, 0.0196;                      %Experiments with 0.0025% equivalents of PhCat
%    0.0025, 4.00, 0.04762;
    0.0025, 2.00, 0.065;
    0.0025, 1.00, 0.272;
    0.0025, 0.75, 0.363;
    0.0025, 0.50, 0.461    
    ];
data005=[                                       %Experiments with 0.005% equivalents of PhCat
%    0.005, 6.00, 0.0392;  
%    0.005, 4.00, 0.08257;
    0.005, 2.00, 0.150;
    0.005, 1.00, 0.33;
    0.005, 0.75, 0.376;
    0.005, 0.50, 0.483
    ];
data02=[                                        %Experiments with 0.02% equivalents of PhCat
%    0.02, 6.00, 0.0180;
%    0.02, 4.00, 0.1105;
    0.02, 2.00, 0.184;
    0.02, 1.00, 0.386;
    0.02, 0.75, 0.4522;
    0.02, 0.50, 0.5845;
    ];
% data05=[                                        %Experiments with 0.05% equivalents of PhCat
%     0.05, 6.00, 0.0566;
%     0.05, 4.00, 0.0698;
%     0.05, 2.00, 0.2308;
%     0.05, 1.00, 0.3691;
%     0.05, 0.75, 0.4083;
%     0.05, 0.50, 0.4853;
%     ];
data07=[                                        %Experiments with 0.07% equivalents of PhCat
%    0.07, 6.00, 0.0610;         
%    0.07, 4.00, 0.0991;
    0.07, 2.00, 0.2752;
    0.07, 1.00, 0.4203;
    0.07, 0.75, 0.4972;
    0.07, 0.50, 0.5908;
    ];    
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Concentration of A
%--------------------------------------------------------------------------
CA0=0.4;                                        %Concentration of A [mol/L]
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
kappa_matrix=160.5;                             %Naperian extintion coefficient of the matrix without photocatalyst [m^-1]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Actinometry results, photon flow rate.
%--------------------------------------------------------------------------
N_Acti=7.7324e-7;                       %Photon flow rate by Actinometry [einstein/s]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Initial guess for LSE method, first column is -log10(molecular diffusitvity); second column phi*k
x0 = [11 0.01];                          %Initial guess for x(1) and x(2)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%End of input data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%Matrix manipulations
%--------------------------------------------------------------------------

datacombined_all=[data0025;data005;data02;data07];

%Fractional conversions
X_exp_all=datacombined_all(:,3)';
X_exp_0025=data0025(:,3)';
X_exp_005=data005(:,3)';
X_exp_02=data02(:,3)';
% X_exp_05=data05(:,3)';
X_exp_07=data07(:,3)';

%Photocatalyst concentration
CP0_exp_all=datacombined_all(:,1).*CA0/100;     
CP0_0025=data0025(:,1).*CA0/100;
CP0_005=data005(:,1).*CA0/100;
CP0_02=data02(:,1).*CA0/100;
% CP0_05=data05(:,1).*CA0/100;
CP0_07=data07(:,1).*CA0/100;

%Flow rate conversion to [m^3 s^-1] and manipulation
Q_exp_all=datacombined_all(:,2)./(10^6*60);     
Q_0025=data0025(:,2)./(10^6*60);                
Q_005=data005(:,2)./(10^6*60);
Q_02=data02(:,2)./(10^6*60);
% Q_05=data05(:,2)./(10^6*60);
Q_07=data07(:,2)./(10^6*60);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Dimensions-conversions and derived geometric variables
%--------------------------------------------------------------------------
Ri_in=Di_in/2;                                  %FEP outer radius in [in]
Ri=Ri_in*0.0254;                                %FEP outer radius in [m]
Ro_in=Do_in/2;                                  %SS internal radius in [in]
Ro=Ro_in*0.0254;                                %SS internal radius in [m]
A_cross=pi*(Ro^2-Ri^2);                         %Cross sectional area in [m^2]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Mean velocity calculations in [m/s]
%--------------------------------------------------------------------------
vmean_data_all=Q_exp_all./A_cross;              
vmean_0025=Q_0025./A_cross;
vmean_005=Q_005./A_cross;
vmean_02=Q_02./A_cross;
% vmean_05=Q_05./A_cross;
vmean_07=Q_07./A_cross;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Naperian attenuation coefficients
%--------------------------------------------------------------------------
kappa_PC_all=alpha_PC.*CP0_exp_all;
kappa_PC_0025=alpha_PC.*CP0_0025;
kappa_PC_005=alpha_PC.*CP0_005;
kappa_PC_02=alpha_PC.*CP0_02;
% kappa_PC_05=alpha_PC.*CP0_05;
kappa_PC_07=alpha_PC.*CP0_07;

kappa_tot_all=kappa_PC_all+kappa_matrix;
kappa_tot_0025=kappa_PC_0025+kappa_matrix;
kappa_tot_005=kappa_PC_005+kappa_matrix;
kappa_tot_02=kappa_PC_02+kappa_matrix;
% kappa_tot_05=kappa_PC_05+kappa_matrix;
kappa_tot_07=kappa_PC_07+kappa_matrix;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Actinometry results, photon flow rate.
%--------------------------------------------------------------------------
N_LDF=N_Acti/(10^-0.07-(10^-1));                    %Photon flow rate by whole LDF [einstein/s]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Calculation of QY and Absorbed photon equivalents
%--------------------------------------------------------------------------
%Experimental data: derived
FA0_exp_all=Q_exp_all.*1000.*CA0;                   %Molar flow rate of A at the inlet [mol/s]
FP0_exp_all=FA0_exp_all.*X_exp_all';                %Molar flow rate of P at the inlet [mol/s]
Quantumyield=FP0_exp_all./N_Acti;                   %External quantum yield
Transmission=exp(-kappa_tot_all*(Ro-Ri));           %Fractional transmission across the annular gap
photon_eq=N_Acti.*(1-Transmission)./FA0_exp_all;    %Absorbed Photon equivalents
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Mesh for PDE
%--------------------------------------------------------------------------
rsteps=100;                                         %number of steps in the radial direction
zsteps=50;                                          %number of steps in the axial direction
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Obtain phi*k and Dm using the LSE method
%--------------------------------------------------------------------------
%Independent data for the model
idata=[kappa_PC_all,kappa_tot_all,vmean_data_all]';

%Declare anonymous function with independent vector 'idata' and dependent
%vector 'x' . x(1) is the -log10(molecular diffusitvity) and x(2) is phi*k
TotX_final=@(x,idata) avg_conversion_out(Ri, Ro, idata(3,:), x(1), x(2),N_LDF, idata(1,:), idata(2,:), L, rsteps, zsteps);

% Define options structure with the desired algorithm
options = optimoptions(@lsqcurvefit, 'Algorithm', 'levenberg-marquardt');
[x,resnorm] = lsqcurvefit(TotX_final, x0, idata, X_exp_all,[],[],options);
% Obtain the parameter values from lsqcurvefit for vector 'x'
pDm= x(1);
phik=x(2);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Calculate predicted values using the obtained parameter values and
% calculate R^2
%--------------------------------------------------------------------------

[TotX_final_prediction] =TotX_final(x, idata);

%Calculate Pearson R^2 (global)
R_sq_all=corrcoef(TotX_final_prediction,X_exp_all);
R_sq_all=R_sq_all(2,1);
R_sq_0025=corrcoef(TotX_final_prediction(1:4),X_exp_0025);
R_sq_0025=R_sq_0025(2,1);
R_sq_005=corrcoef(TotX_final_prediction(5:8),X_exp_005);
R_sq_005=R_sq_005(2,1);
R_sq_020=corrcoef(TotX_final_prediction(9:12),X_exp_02);
R_sq_020=R_sq_020(2,1);
R_sq_070=corrcoef(TotX_final_prediction(13:16),X_exp_07);
R_sq_070=R_sq_070(2,1);
straightline=[0 1];                         %For correlation R^2=1
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Figures
%--------------------------------------------------------------------------
figure(1); %Plot experimental data; X(Q)
plot(data0025(:,2), data0025(:,3), 'o', 'MarkerFaceColor', "#0072BD");
hold on;
plot(data005(:,2), data005(:,3), 's', 'MarkerFaceColor', "#D95319");
plot(data02(:,2), data02(:,3), 'd', 'MarkerFaceColor', "#EDB120");
% plot(data05(:,2), data05(:,3), 'v', 'MarkerFaceColor', "#77AC30");
plot(data07(:,2), data07(:,3), '^', 'MarkerFaceColor', "#7E2F8E");
xlabel('$Q$ [mL min$^{-1}$]', 'Interpreter', 'latex');
ylabel('$\left\langle X_A \right\rangle$ [-]', 'Interpreter', 'latex');
legend('0.0025%', '0.0050%', '0.0200%','0.0700%');

figure(2); %Plot experimental vs. prediction (Parity plot) 
plot(TotX_final_prediction(1:4), X_exp_0025,'o', 'MarkerFaceColor', "#0072BD");         %Experimental set 1
hold on;
plot(TotX_final_prediction(5:8), X_exp_005, 's', 'MarkerFaceColor', "#D95319");    %Experimental set 2
plot(TotX_final_prediction(9:12), X_exp_02, 'd', 'MarkerFaceColor', "#EDB120");        %Experimental set 3
% plot(TotX_final_prediction(19:24), X_exp_05, 'v', 'MarkerFaceColor', "#77AC30");        %Experimental set 4
plot(TotX_final_prediction(13:16), X_exp_07,'^', 'MarkerFaceColor', "#7E2F8E");         %Experimental set 5
plot(straightline,straightline,'color', 'black');
xlabel('PDE Model $\left\langle X_A \right\rangle$ [-]', 'Interpreter', 'latex');
ylabel('Experimental $\left\langle X_A \right\rangle$ [-]', 'Interpreter', 'latex');
legend('0.0025%','0.0050%', '0.0200%','0.0700%');

% Calculate predicted values for smooth plot X(Q)
smooth_steps=30;                    %number of points oto evaluate per data set
%Entend line to 30% over and 40% under the experimental range
vmean_data_dummy=linspace(vmean_data_all(1)*1.3, vmean_data_all(4)*0.6, smooth_steps);
%Flow rate for whole range
Q_mlmin_dummy=A_cross*10^6*60.*vmean_data_dummy;

kappa_PC_0025_dummy=repmat(kappa_PC_0025(1),smooth_steps,1)';
kappa_tot_0025_dummy=repmat(kappa_tot_0025(1),smooth_steps,1)';
dummy_idata_0025=[kappa_PC_0025_dummy;kappa_tot_0025_dummy;vmean_data_dummy];

kappa_PC_005_dummy=repmat(kappa_PC_005(1),smooth_steps,1)';
kappa_tot_005_dummy=repmat(kappa_tot_005(1),smooth_steps,1)';
dummy_idata_005=[kappa_PC_005_dummy;kappa_tot_005_dummy;vmean_data_dummy];

kappa_PC_02_dummy=repmat(kappa_PC_02(1),smooth_steps,1)';
kappa_tot_02_dummy=repmat(kappa_tot_02(1),smooth_steps,1)';
dummy_idata_02=[kappa_PC_02_dummy;kappa_tot_02_dummy;vmean_data_dummy];

% kappa_PC_05_dummy=repmat(kappa_PC_05(1),smooth_steps,1)';
% kappa_tot_05_dummy=repmat(kappa_tot_05(1),smooth_steps,1)';
% dummy_idata_05=[kappa_PC_05_dummy;kappa_tot_05_dummy;vmean_data_dummy];

kappa_PC_07_dummy=repmat(kappa_PC_07(1),smooth_steps,1)';
kappa_tot_07_dummy=repmat(kappa_tot_07(1),smooth_steps,1)';
dummy_idata_07=[kappa_PC_07_dummy;kappa_tot_07_dummy;vmean_data_dummy];

%Evaluate function using x obtained by least square method
[TotX_final_prediction_smooth_0025] =TotX_final(x, dummy_idata_0025);
[TotX_final_prediction_smooth_005] =TotX_final(x, dummy_idata_005);
[TotX_final_prediction_smooth_02] =TotX_final(x, dummy_idata_02);
% [TotX_final_prediction_smooth_05] =TotX_final(x, dummy_idata_05);
[TotX_final_prediction_smooth_07] =TotX_final(x, dummy_idata_07);

%Plot experimental data with prediction curve
figure(3);
plot(data0025(:,2), data0025(:,3), 'o', 'MarkerFaceColor', "#0072BD");
hold on;
plot(data005(:,2), data005(:,3), 's', 'MarkerFaceColor', "#D95319");
plot(data02(:,2), data02(:,3), 'd', 'MarkerFaceColor', "#EDB120");
% plot(data05(:,2), data05(:,3), 'v', 'MarkerFaceColor', "#77AC30");
plot(data07(:,2), data07(:,3), '^', 'MarkerFaceColor', "#7E2F8E");
plot(Q_mlmin_dummy, TotX_final_prediction_smooth_0025, 'color', "#0072BD");
plot(Q_mlmin_dummy, TotX_final_prediction_smooth_005, 'color', "#D95319");
plot(Q_mlmin_dummy, TotX_final_prediction_smooth_02, 'color', "#EDB120");
% plot(Q_mlmin_dummy, TotX_final_prediction_smooth_05, 'color', "#77AC30");
plot(Q_mlmin_dummy, TotX_final_prediction_smooth_07, 'color', "#7E2F8E");
xlabel('$Q$ [mL min$^{-1}$]', 'Interpreter', 'latex');
ylabel('$$\left\langle X_A \right\rangle$$ [-]', 'Interpreter', 'latex');
legend('0.0025%', '0.0050%', '0.0200%', '0.0700%');

%Plot experimental data in terms of photon equivalents
figure(4);
plot(photon_eq(1:4), X_exp_0025, 'o', 'MarkerFaceColor', "#0072BD");
hold on;
plot(photon_eq(5:8), X_exp_005, 's', 'MarkerFaceColor', "#D95319");
plot(photon_eq(9:12), X_exp_02, 'd', 'MarkerFaceColor', "#EDB120");
% plot(photon_eq(19:24), X_exp_05, 'v', 'MarkerFaceColor', "#77AC30");
plot(photon_eq(13:16), X_exp_07, '^', 'MarkerFaceColor', "#7E2F8E");
xlabel('$\eta_{eq,abs}$ [einstein mol$^{-1}$]', 'Interpreter', 'latex');
ylabel('$$\left\langle X_A \right\rangle$$ [-]', 'Interpreter', 'latex');
legend('0.0025%','0.0050%', '0.0200%', '0.0700%');

% Save the constants to a .mat file
save('constants.mat', 'pDm', 'phik');