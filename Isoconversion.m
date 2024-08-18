%--------------------------------------------------------------------------
% Title: CAP-Flow system PDE Model: Isoconversion curves
% Author: Sergio Carrillo De Hert
% Affiliation: University College Dublin
% Last modified: March 05, 2024
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Description: The aim of to create the tau(C_PC) iso-conversion curves for
% the PDE and plu-flow models
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
%           1. Iso-conversion curves to be obtained
%           2. C_PC range to be tested
%           3. Geometry of the CAP-Flow system
%           4. Photon absorption properties of matrix and PC
%           5. Photon flow rate (by actinometry)
%           6. Kinetic contants and molecular diffusion coefficient (from
%           Fitter)
% The output includes:
%           1. Plot of X(Dm) for varying Q for PDE and plug-flow models
%--------------------------------------------------------------------------

clc;
clear;

%--------------------------------------------------------------------------
%Isoconversion curves
%--------------------------------------------------------------------------
X_CAPFLOW_isocurves = [0.9, 0.8, 0.6, 0.4];  
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Concentration of A  and photocatalyst and C_PC space
%--------------------------------------------------------------------------
CA0=0.4;                                            %Concentration of A [mol/L]
noisopoints=45;                                     %number of PhCat concentrations per isoconversion-curve
CP0_eq=linspace(0.0025,0.1,noisopoints);            %create linear space of PC equivalents in percentage
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

% Load the constants from the .mat file
load('constants.mat', 'pDm', 'phik');
% pDm=11;                                             %Molecular diffusion coefficient in [m^2/s]
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
V_R=A_cross*L;                                      %Volume of reactor in [m^3]
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Naperian extintion coefficient
%--------------------------------------------------------------------------
CP0=CP0_eq.*CA0./100;                               %Conc. PhCat in [mol/L]
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
% Repeat each element of X_CAPFLOW 20 times in each column
%--------------------------------------------------------------------------
for i = 1:numel(X_CAPFLOW_isocurves)
    X_outlet(:, i) = repmat(X_CAPFLOW_isocurves(i), noisopoints, 1);
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%PDE solution
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Preallocate variables for speed
%--------------------------------------------------------------------------
vmean_solution=zeros(noisopoints,length(X_CAPFLOW_isocurves)); 
Q=zeros(noisopoints,length(X_CAPFLOW_isocurves)); 
t_res=zeros(noisopoints,length(X_CAPFLOW_isocurves)); 
t_res_hr=zeros(noisopoints,length(X_CAPFLOW_isocurves)); 
Productivity=zeros(noisopoints,length(X_CAPFLOW_isocurves)); 
%--------------------------------------------------------------------------

for j=1:length(X_CAPFLOW_isocurves)
    for i = 1:length(CP0_eq)
        if i==1
            % Define bounds based on X_CAPFLOW_isocurves
            if X_CAPFLOW_isocurves(j) >=0.8
                x0 = 3.29*10^-5;
            else 
                x0 = 3.5*10^-4;
            end
        else
            x0=vmean_solution(i-1,j);
        end

        % Declare anonymous function for vmean_data
        TotX_final = @(vmean_data) avg_conversion_out(Ri, Ro, vmean_data, pDm, phik, N_LDF, kappa_PC(i), kappa_tot(i), L, rsteps, zsteps);
        % Solve for vmean_data
        vmean_solution(i,j) = fzero(@(vmean_data) TotX_final(vmean_data) - X_outlet(i,j), x0);
        Q(i,j)=vmean_solution(i,j).*A_cross;
        t_res(i,j)=V_R/Q(i,j);
        t_res_hr(i,j)=t_res(i,j)./3600;
        Productivity(i,j)=Q(i,j).*1000*CA0.*X_outlet(i,j)*10^6;
    end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Ideal plug flow
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Preallocate variables for speed
%--------------------------------------------------------------------------
Q_plug=zeros(noisopoints,length(X_CAPFLOW_isocurves)); 
t_res_plug=zeros(noisopoints,length(X_CAPFLOW_isocurves)); 
t_res_hr_plug=zeros(noisopoints,length(X_CAPFLOW_isocurves)); 
Productivity_plug=zeros(noisopoints,length(X_CAPFLOW_isocurves)); 
%--------------------------------------------------------------------------

for j=1:length(X_CAPFLOW_isocurves)
    for i = 1:length(CP0_eq)
        Q_plug(i,j) = (phik/log(1-X_outlet(i,j)))*N_LDF*(1/10)*(10^L-1)*(exp(-kappa_tot(i)*(Ro-Ri))-1);
        t_res_plug(i,j)=V_R/Q_plug(i,j);
        t_res_hr_plug(i,j)=t_res_plug(i,j)./3600;
        Productivity_plug(i,j)=Q_plug(i,j).*1000*CA0.*X_outlet(i,j)*10^6;
    end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Figures
%--------------------------------------------------------------------------

% Define a color map for the curves
colors = lines(length(X_CAPFLOW_isocurves));
figure(1);
hold on;
% Get handle to current axes
ax = gca;

% Turn on the box around the plot
ax.Box = 'on';
% Plot continuous lines with legends
for i = 1:length(X_CAPFLOW_isocurves)
    color = colors(i, :);
    % Plot t_res_hr data with continuous lines
    h1=plot(CP0_eq/100, t_res_hr(:,i), '-', 'Color', color);
    h2=plot(CP0_eq/100, t_res_hr_plug(:,i), '--', 'Color', color);
    % Add legend label only for the continuous line
    legendLabels{i} = sprintf('%0.2f', X_CAPFLOW_isocurves(i));
    % Exclude the dashed line from the legend
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
hold off;
% Create legend for continuous lines
ylabel('$\tau$ [h]', 'Interpreter', 'latex');
xlabel('$C_{PC}$ [equiv.]', 'Interpreter', 'latex');
lgd=legend(legendLabels, 'Interpreter', 'latex');
% Add a title to the legend
title(lgd, '$\left\langle X_A \right\rangle$ [-]', 'Interpreter', 'latex');
hold off;

figure(2);
hold on;
for i=1:length(X_CAPFLOW_isocurves)
    color = colors(i, :);
    h3=plot(CP0_eq, Productivity(:,i), '-', 'Color', color);
    h4=plot(CP0_eq, Productivity_plug(:,i), '--', 'Color', color);
    % Add legend label only for the continuous line
    legendLabels2{i} = sprintf('%0.2f', X_CAPFLOW_isocurves(i));
    % Exclude the dashed line from the legend
    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
ylabel('$F_P$ [$\mu$mol/s]', 'Interpreter', 'latex');
xlabel('$C_{PC}$ [\%]', 'Interpreter', 'latex');
legend(legendLabels2, 'Interpreter', 'latex');