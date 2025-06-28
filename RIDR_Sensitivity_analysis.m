%% RIDR Sensitivity Analysis for Demolition Loss Assessment
% =========================================================================
% PROBABILISTIC DEMOLITION LOSS ASSESSMENT - SENSITIVITY TO RIDR PARAMETERS
% =========================================================================
%
% Description:
%   This script performs sensitivity analysis on Residual Inter-story Drift
%   Ratio (RIDR) parameters for demolition loss%   assessment. 
%   The analysis evaluates how variations in demolition threshold
%   (theta) and uncertainty (beta) affect Expected Annual Loss (EAL) due to
%   building demolition following earthquake damage.
%
% Author: Shiva Baddipalli
% Email: shivalinga.baddipalli@usu.edu
% Last Updated: June 27, 2025
%
% Corresponding Author:
% Dr. Mohsen Zaker Esteghamati (Assistant Professor): mohsen.zaker@usu.edu
% Department of Civil and Environmental Engineering
% Utah State University, Logan, UT, United States
%
% Reference:
%   Ramirez, C.M., and Miranda, E. (2012). Significance of residual drifts
%   in building earthquake loss estimation. Earthquake Engineering and
%   Structural Dynamics, 41(11), 1477-1493.
%
% Input Files Required:
%   - GuanDataBase-IMs.csv: Intensity measures (Sa) for 240 ground motions
%   - Building geometry files in Guan_database/BuildingDesigns/
%   - Residual story drift files (RSDR) in Guan_database/StructuralResponses/
%   - USGSHazard_data.csv: Seismic hazard curves for different periods
%
% Output:
%   - Expected Annual Loss (EAL) for demolition [$]
%   - Normalized EAL for demolition (as fraction of replacement cost) [-]
%   - Sensitivity grid results for theta and beta parameter combinations
%
% Usage:
%   Run this script directly. Results are saved in Excel files within the
%   "RIDR_Sensitivity_Grid_Data" folder for each parameter combination.
%
% =========================================================================

%% Initialize Workspace
% Clear workspace and command window for clean execution
clc; 
clearvars; 
close all;
clear;

%% Data Input and Configuration
% -------------------------------------------------------------------------
% Load ground motion intensity measures and building information
% -------------------------------------------------------------------------

% File path for intensity measure database
IM_FilePath = 'GuanDataBase-IMs.csv';

% Read intensity measure data
fprintf('Loading ground motion intensity measures...\n');
Data = readmatrix(IM_FilePath);

% Extract relevant data columns
IM_GM = Data(:, 3:242);              % IM values for all 240 ground motions (building-wise)
ID = Data(:, 1);                     % Building ID numbers
IM_GM_trans = IM_GM';                % Transposed IM matrix for easier processing
Target_time_period = Data(:, 2);     % Fundamental periods of buildings [seconds]

% Define spectral acceleration range for integration
SaCalc = [0.001, 0.01, 0.05, 0.1:0.01:6]';  % Sa values [g]

fprintf('Loaded data for %d buildings with %d ground motions each.\n', length(ID), size(IM_GM,2));

%% RIDR Sensitivity Parameter Definition
% -------------------------------------------------------------------------
% Define parameter ranges for sensitivity analysis
% -------------------------------------------------------------------------

% Demolition threshold parameters (median RIDR causing demolition)
theta_hat_Demolish = [0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095];  % [rad]

% Demolition uncertainty parameters (logarithmic standard deviation)
beta_hat_Demolish = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];           % [-]

fprintf('Sensitivity analysis parameters:\n');
fprintf('  - Theta (demolition threshold): %d values from %.3f to %.3f\n', ...
        length(theta_hat_Demolish), min(theta_hat_Demolish), max(theta_hat_Demolish));
fprintf('  - Beta (uncertainty): %d values from %.1f to %.1f\n', ...
        length(beta_hat_Demolish), min(beta_hat_Demolish), max(beta_hat_Demolish));
fprintf('  - Total parameter combinations: %d\n\n', ...
        length(theta_hat_Demolish) * length(beta_hat_Demolish));

%% Initialize Output Variables
% -------------------------------------------------------------------------
% Pre-allocate arrays for computational efficiency
% -------------------------------------------------------------------------

dlambdaSa = zeros(length(SaCalc), length(ID));     % Hazard density function
EAL_Demo = zeros(length(ID), 1);                   % Expected Annual Loss - Demolition [$]
Norm_EAL_Demo = zeros(length(ID), 1);              % Normalized EAL - Demolition [-]

%% Main Sensitivity Analysis Loop
% =========================================================================
% Nested loops for parameter sensitivity analysis
% =========================================================================

fprintf('Starting sensitivity analysis...\n');
total_combinations = length(beta_hat_Demolish) * length(theta_hat_Demolish);
current_combination = 0;

% Loop over beta values (demolition uncertainty)
for b_dem = 1:length(beta_hat_Demolish)
    
    % Loop over theta values (demolition threshold)
    for t_dem = 1:length(theta_hat_Demolish)
        
        current_combination = current_combination + 1;
        fprintf('\nProcessing combination %d/%d: Theta=%.3f, Beta=%.1f\n', ...
                current_combination, total_combinations, ...
                theta_hat_Demolish(t_dem), beta_hat_Demolish(b_dem));
        
        % Reset EAL arrays for each parameter combination
        EAL_Demo = zeros(length(ID), 1);
        Norm_EAL_Demo = zeros(length(ID), 1);
        
        %% Building-wise Analysis Loop
        % -----------------------------------------------------------------
        % Process each building individually
        % -----------------------------------------------------------------
        
        for bd = 1:length(ID)
            
            %% File Path Definition and Data Loading
            % -------------------------------------------------------------
            % Construct file paths for current building
            % -------------------------------------------------------------
            
            % Residual story drift ratio file path
            fileRSDR = strcat(['Guan_database\Structural', ...
                'Responses\EDPsUnder240GMs\ResidualStoryDrift'], '\', ...
                'Building_', num2str(ID(bd),'%.1d'), '_RSDR.csv');
            
            % Building geometry file path
            fileGeometry = strcat('Guan_database\BuildingDesigns', '\', ...
                'Building_', num2str(ID(bd),'%.1d'), '\', 'Geometry.csv');
            
            % Load structural response data
            RSDR = readmatrix(fileRSDR);                % Residual story drift ratios
            Geometry = readmatrix(fileGeometry);         % Building geometry parameters
            
            %% Process Residual Story Drift Data
            % -------------------------------------------------------------
            % Extract maximum RSDR for each ground motion
            % -------------------------------------------------------------
            
            Max_RSDR = zeros(1, length(RSDR));
            for gm = 1:length(RSDR)
                Max_RSDR(1, gm) = max(abs(RSDR(:, gm)));  % Maximum absolute RSDR
            end
            
            %% Building Geometry and Cost Parameters
            % -------------------------------------------------------------
            % Extract building dimensions and calculate replacement cost
            % -------------------------------------------------------------
            
            Num_stories = Geometry(1, 1);               % Number of stories
            Bay_width = Geometry(1, 6);                 % Bay width [ft]
            Num_bays = 5;                               % Number of bays (assumed)
            
            % Calculate floor area and total replacement cost
            Floor_area = (Bay_width * Num_bays) * (Bay_width * Num_bays);  % [ft²]
            Per_sqft_Cost = 250;                        % Cost per square foot [$/ft²]
            Tot_Replac_Cost = Num_stories * Floor_area * Per_sqft_Cost;    % Total replacement cost [$]
            
            % Demolition cost components
            DemoLoss = 0.1 * Tot_Replac_Cost;          % Demolition cost (10% of replacement)
            Tot_Loss_Demo = Tot_Replac_Cost + DemoLoss; % Total loss if demolished
            
            %% Demolition Probability Assessment
            % -------------------------------------------------------------
            % Calculate probability of demolition based on RIDR
            % -------------------------------------------------------------
            
            % Define RSDR range for probability density function
            x_RSDR_pdf = (0.0002:0.0002:0.06)';         % RSDR values for integration
            
            % Demolition fragility curve (lognormal CDF)
            % Based on Ramirez and Miranda (2012)
            PD_RSDR = normcdf((log(x_RSDR_pdf / theta_hat_Demolish(t_dem))) / ...
                             beta_hat_Demolish(b_dem));
            
            % Initialize demolition probability array
            PD_IM_NC = zeros(length(SaCalc), 1);
            
            %% Demand Model Development and Integration
            % -------------------------------------------------------------
            % Develop RSDR demand model and calculate demolition probability
            % -------------------------------------------------------------
            
            % Fit log-linear demand model: ln(RSDR) = a + b*ln(Sa) + error
            md_RSDR = fitlm(log(IM_GM(bd, :)), log(Max_RSDR(1, :)));
            
            % Extract demand model parameters
            RSDR_DemandPar(1:2) = md_RSDR.Coefficients{1:2, 1};  % [intercept, slope]
            RSDR_DemandPar(3) = md_RSDR.RMSE;                     % Model uncertainty
            SIGMA_RSDR = RSDR_DemandPar(3);
            
            % Calculate demolition probability for each Sa level
            for im = 1:length(SaCalc)
                
                % Mean of log(RSDR) given Sa
                MD_RSDR = RSDR_DemandPar(2) * log(SaCalc(im)) + RSDR_DemandPar(1);
                
                % Probability density function of RSDR given Sa
                dP_RIDR_IM = lognpdf(x_RSDR_pdf, MD_RSDR, SIGMA_RSDR);
                
                % Convolution of demolition fragility with RSDR distribution
                PD_RSDR_x_dP_RSDR_IM = PD_RSDR(:, 1) .* dP_RIDR_IM(:, 1);
                
                % Integrate to get demolition probability
                PD_IM_NC(im, 1) = trapz(x_RSDR_pdf, PD_RSDR_x_dP_RSDR_IM);
            end
            
            PD_IM_NC = PD_IM_NC';  % Transpose for consistency
            
            %% Demolition Loss Calculation
            % -------------------------------------------------------------
            % Calculate expected loss due to demolition
            % -------------------------------------------------------------
            
            % Loss given demolition occurs
            Loss_Demo(1, :) = Tot_Loss_Demo * PD_IM_NC(1, :);
            
            % Normalized loss (fraction of replacement cost)
            Norm_Loss_Demo(1, :) = Loss_Demo(1, :) ./ Tot_Replac_Cost;
            
            %% Seismic Hazard Processing
            % -------------------------------------------------------------
            % Load and process seismic hazard data for current building
            % -------------------------------------------------------------
            
            % Load USGS hazard data
            Hazard_file = 'USGSHazard_data.csv';
            HazardData = readmatrix(Hazard_file);
            AllLamda_IM = HazardData(:, 2:12);          % Hazard rates for different periods
            IM_values = HazardData(:, 1);               % Sa values for hazard curve
            
            % Extract hazard rates for different periods
            hazard_columns = num2cell(AllLamda_IM, 1);
            [a, b, c, d, e, f, g, h, i, j] = hazard_columns{:};
            
            % Interpolate hazard curve for building's fundamental period
            Lamda_IM_target = zeros(length(IM_values), length(Target_time_period));
            T_target = Target_time_period(bd, 1);
            
            % Period-based interpolation of hazard curves
            if T_target >= 0.1 && T_target < 0.2
                for nn = 1:length(IM_values)
                    Lamda_IM_target(nn, bd) = interp1([0.1, 0.2], ...
                        [log(a(nn, 1)), log(b(nn, 1))], T_target, 'linear');
                end
            elseif T_target >= 0.2 && T_target < 0.3
                for nn = 1:length(IM_values)
                    Lamda_IM_target(nn, bd) = interp1([0.2, 0.3], ...
                        [log(b(nn, 1)), log(c(nn, 1))], T_target, 'linear');
                end
            elseif T_target >= 0.3 && T_target < 0.5
                for nn = 1:length(IM_values)
                    Lamda_IM_target(nn, bd) = interp1([0.3, 0.5], ...
                        [log(c(nn, 1)), log(d(nn, 1))], T_target, 'linear');
                end
            elseif T_target >= 0.5 && T_target < 0.75
                for nn = 1:length(IM_values)
                    Lamda_IM_target(nn, bd) = interp1([0.5, 0.75], ...
                        [log(d(nn, 1)), log(e(nn, 1))], T_target, 'linear');
                end
            elseif T_target >= 0.75 && T_target < 1
                for nn = 1:length(IM_values)
                    Lamda_IM_target(nn, bd) = interp1([0.75, 1], ...
                        [log(e(nn, 1)), log(f(nn, 1))], T_target, 'linear');
                end
            elseif T_target >= 1 && T_target < 2
                for nn = 1:length(IM_values)
                    Lamda_IM_target(nn, bd) = interp1([1, 2], ...
                        [log(f(nn, 1)), log(g(nn, 1))], T_target, 'linear');
                end
            elseif T_target >= 2 && T_target < 3
                for nn = 1:length(IM_values)
                    Lamda_IM_target(nn, bd) = interp1([2, 3], ...
                        [log(g(nn, 1)), log(h(nn, 1))], T_target, 'linear');
                end
            elseif T_target >= 3 && T_target < 4
                for nn = 1:length(IM_values)
                    Lamda_IM_target(nn, bd) = interp1([3, 4], ...
                        [log(h(nn, 1)), log(i(nn, 1))], T_target, 'linear');
                end
            elseif T_target >= 4 && T_target < 5
                for nn = 1:length(IM_values)
                    Lamda_IM_target(nn, bd) = interp1([4, 5], ...
                        [log(i(nn, 1)), log(j(nn, 1))], T_target, 'linear');
                end
            end
            
            %% Hazard Curve Processing and Differentiation
            % -------------------------------------------------------------
            % Convert hazard curve to probability density function
            % -------------------------------------------------------------
            
            % Create hazard curve for current building
            hazardCurveSa = [IM_values, exp(Lamda_IM_target(:, bd))];
            
            % Numerical differentiation to get hazard density
            PDF_temp = zeros(length(SaCalc), 1);
            hh = 0.0000001;  % Small increment for numerical differentiation
            
            for kk = 1:length(SaCalc)
                % Two-point finite difference formula
                fxp = interp1(hazardCurveSa(:, 1), hazardCurveSa(:, 2), ...
                             SaCalc(kk) + hh, 'spline');
                fxn = interp1(hazardCurveSa(:, 1), hazardCurveSa(:, 2), ...
                             SaCalc(kk) - hh, 'spline');
                PDF_temp(kk) = abs((fxp - fxn)) / (2 * hh);
            end
            
            % Smooth the hazard density function
            dlambdaSa(:, bd) = smooth(PDF_temp);
            
            %% Expected Annual Loss (EAL) Computation
            % -------------------------------------------------------------
            % Calculate EAL by integrating loss with hazard density
            % -------------------------------------------------------------
            
            % Convolution of loss with hazard density
            Conv_loss_Demo = Loss_Demo(1, :)' .* dlambdaSa(:, bd);
            Conv_Norm_loss_Demo = Norm_Loss_Demo(1, :)' .* dlambdaSa(:, bd);
            
            % Integration over Sa range to get EAL
            EAL_Demo(bd, 1) = trapz(SaCalc, Conv_loss_Demo);
            Norm_EAL_Demo(bd, 1) = trapz(SaCalc, Conv_Norm_loss_Demo);
            
        end  % End building loop
        
        %% Results Compilation and Output
        % -----------------------------------------------------------------
        % Prepare and save results for current parameter combination
        % -----------------------------------------------------------------
        
        % Compile results matrix
        All_Output = zeros(length(ID), 3);
        All_Output(:, 1) = ID(:, 1);                    % Building ID
        All_Output(:, 2) = EAL_Demo(:, 1);              % EAL - Demolition [$]
        All_Output(:, 3) = Norm_EAL_Demo(:, 1);         % Normalized EAL - Demolition [-]
        
        % Generate output filename with parameter values
        excelFileName = sprintf('RIDR_Sensitivity_Grid_Data\\RIDR_Sensitivity_DemoLoss_Theta_%.3f_Beta_%.1f.xlsx', ...
                               theta_hat_Demolish(t_dem), beta_hat_Demolish(b_dem));
        
        % Save results to Excel file (UNCOMMENT THE LINES BELOW TO ENABLE FILE SAVING)
        % xlswrite(excelFileName, All_Output, 'Sheet1');
        % fprintf('  → Results saved: %s\n', excelFileName);
        
    end  % End theta loop
    
end  % End beta loop


%% End of Script
