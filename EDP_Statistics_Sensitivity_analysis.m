%% EAL Assessment Using Weighted EDP Statistics with Regression-Based Alpha
% =========================================================================
% SEISMIC LOSS ASSESSMENT USING PERIOD-DEPENDENT ALPHA REGRESSION MODEL
% =========================================================================
%
% Description:
%   This script implements seismic loss assessment using weighted EDP statistics
%   where the weighting factor (alpha) is determined from a regression model
%   based on building fundamental period. The approach applies optimized alpha
%   values from previous calibration to compute normalized Expected Annual Loss
%   (EAL) using HAZUS assembly-based methodology.
%
% Methodology:
%   - Uses regression equation: log(α) = f(log(T)) where T is fundamental period
%   - Applies weighted EDP statistics: Weighted_EDP = α × P65 + (1-α) × P35
%   - Calculates component-level losses (structural, non-structural)
%   - Includes demolition probability based on residual drift
%   - Integrates with seismic hazard to compute EAL values
%
% Author: Shiva Baddipalli
% Institution: Utah State University, USA
% Email: shivalinga.baddipalli@usu.edu
% Last Updated: June 27, 2025
%
% Reference:
%   FEMA (2003). Multi-hazard Loss Estimation Methodology - Earthquake Model
%   HAZUS®MH MR4 Technical Manual. Federal Emergency Management Agency.
%
% Input Files Required:
%   - GuanDataBase-IMs.csv: Intensity measures for 240 ground motions
%   - Building_Info_URGENT.xlsx: Building information including periods
%   - Building response files (PSDR, PFA, RSDR) in Guan_database/
%   - Fragility and cost data files in HAZUS/
%   - USGSHazard_data.csv: Seismic hazard curves
%
% Output:
%   - HAZUS_NormEAL_Weighted_Testing_35_65_10-07-24.xlsx: Final EAL results
%   - Console output showing execution time and completion status
%
% Usage:
%   Run this script directly. Ensure all required data files are in the
%   correct directory structure. The regression model determines alpha
%   automatically for each building based on its fundamental period.
%
% =========================================================================

%% Initialize Workspace and Timer
% Clear workspace and start timing the analysis
clc;                                     % Clear command window
clearvars;                               % Clear all variables from workspace
close all;                               % Close all figure windows
clear;                                   % Clear workspace (redundant but ensures clean start)
tic;                                     % Start timing the analysis

fprintf('=================================================================\n');
fprintf('EAL ASSESSMENT USING WEIGHTED EDP STATISTICS\n');
fprintf('REGRESSION-BASED ALPHA MODEL (25TH-75TH PERCENTILES)\n');
fprintf('=================================================================\n\n');

%% Data Input and Configuration
% -------------------------------------------------------------------------
% Load input data and initialize analysis parameters
% -------------------------------------------------------------------------

% Load ground motion intensity measures from Guan database
IM_FilePath = 'GuanDataBase-IMs.csv';
fprintf('Loading ground motion intensity measures...\n');
Data = readmatrix(IM_FilePath);          % Read entire CSV file into matrix

% Extract data components from loaded matrix
IM_GM = Data(:, 3:242);                  % IM values for all 240 ground motions (columns 3-242, building-wise)
ID = Data(:, 1);                         % Building ID numbers (column 1)
Target_time_period = Data(:, 2);         % Fundamental periods of buildings in seconds (column 2)
SaCalc = [0.001, 0.01, 0.05, 0.1:0.01:6]';  % Spectral acceleration range for integration [g]

% Display loaded data summary for user verification
fprintf('Data loaded successfully:\n');
fprintf('  - Buildings: %d\n', length(ID));
fprintf('  - Ground motions per building: %d\n', size(IM_GM, 2));
fprintf('  - Sa integration points: %d\n\n', length(SaCalc));

%% Initialize Output Variables
% -------------------------------------------------------------------------
% Pre-allocate arrays for computational efficiency
% -------------------------------------------------------------------------

% Hazard-related arrays (will be computed for each building)
dlambdaSa = zeros(length(SaCalc), length(ID));     % Hazard density function matrix

% Expected Annual Loss arrays (absolute dollar values)
EAL_Tot = zeros(length(ID), 1);                    % Total EAL [$]
EAL_Demo = zeros(length(ID), 1);                   % Demolition EAL [$]
EAL_Struc = zeros(length(ID), 1);                  % Structural EAL [$]
EAL_NonStrucDrift = zeros(length(ID), 1);          % Non-structural drift-sensitive EAL [$]
EAL_NonStrucAcc = zeros(length(ID), 1);            % Non-structural acceleration-sensitive EAL [$]

% Normalized Expected Annual Loss arrays (as fraction of replacement cost)
Norm_EAL_Tot = zeros(length(ID), 1);               % Total normalized EAL [-]
Norm_EAL_Demo = zeros(length(ID), 1);              % Demolition normalized EAL [-]
Norm_EAL_Struc = zeros(length(ID), 1);             % Structural normalized EAL [-]
Norm_EAL_NonStrucDrift = zeros(length(ID), 1);     % Non-structural drift-sensitive normalized EAL [-]
Norm_EAL_NonStrucAcc = zeros(length(ID), 1);       % Non-structural acceleration-sensitive normalized EAL [-]

% Alpha values storage (computed from regression model)
alpha_values = zeros(length(ID), 1);               % Alpha values for each building

% Demolition parameters (fixed values from Ramirez and Miranda 2012)
theta_hat_Demolish = 0.015;                       % Median RIDR threshold for demolition decision [rad]
beta_hat_Demolish = 0.3;                          % Logarithmic standard deviation for demolition fragility [-]

fprintf('Analysis parameters:\n');
fprintf('  - Weighted EDP approach: 35th-65th percentiles\n');
fprintf('  - Alpha determination: Regression model based on period\n');
fprintf('  - Demolition threshold (theta): %.3f rad\n', theta_hat_Demolish);
fprintf('  - Demolition uncertainty (beta): %.1f\n\n', beta_hat_Demolish);

%% Main Building Loop - Loss Assessment
% =========================================================================
% Iterate through each building to calculate EAL using weighted EDP approach
% =========================================================================

fprintf('Starting building-wise loss assessment...\n\n');

for bd = 1:length(ID)
    
    % Display current progress to user
    fprintf('Processing Building %d/%d (ID: %d)...\n', bd, length(ID), ID(bd));
    
    %% File Path Definition and Data Loading
    % ---------------------------------------------------------------------
    % Load structural response and geometry data for current building
    % ---------------------------------------------------------------------
    
    % Load building information file (contains fundamental periods and other parameters)
    csvFilePath1 = readmatrix('Building_Info.xlsx');
    
    % Construct file paths for current building's response data
    % Note: File naming convention follows Guan database structure
    filePSDR = strcat(['Guan_database\Structura', ...
        'lResponses\EDPsUnder240GMs\PeakStoryDrift'], '\', ...
        'Building_', num2str(ID(bd), '%.1d'), '_PSDR.csv');
    filePFA = strcat(['Guan_database\Structural', ...
        'Responses\EDPsUnder240GMs\PeakFloorAcceleration'], '\', ...
        'Building_', num2str(ID(bd), '%.1d'), '_PFA.csv');
    fileRSDR = strcat(['Guan_database\Structural', ...
        'Responses\EDPsUnder240GMs\ResidualStoryDrift'], '\', ...
        'Building_', num2str(ID(bd), '%.1d'), '_RSDR.csv');
    fileGeometry = strcat('Guan_database\BuildingDesigns', '\', ...
        'Building_', num2str(ID(bd), '%.1d'), '\', 'Geometry.csv');
    
    % Load structural response data matrices
    PSDR = readmatrix(filePSDR);                   % Peak story drift ratios [rad] - matrix: floors × ground motions
    PFA = readmatrix(filePFA);                     % Peak floor accelerations [in/s²] - matrix: floors × ground motions  
    RSDR = readmatrix(fileRSDR);                   % Residual story drift ratios [rad] - matrix: floors × ground motions
    Geometry = readmatrix(fileGeometry);           % Building geometry parameters
    
    %% Process Engineering Demand Parameters (EDPs)
    % ---------------------------------------------------------------------
    % Calculate statistical measures for each EDP type across all ground motions
    % ---------------------------------------------------------------------
    
    % Convert PFA from inches/s² to acceleration in g units
    PFA = PFA * 0.0254 * 0.10197;                 % Conversion factor: in/s² to g
    
    % Get number of ground motions for current building
    num_gms = size(PSDR, 2);                      % Number of columns = number of ground motions
    
    % Initialize statistical arrays for Peak Story Drift Ratio (PSDR)
    Max_PSDR = zeros(1, num_gms);                 % Maximum PSDR across all floors for each GM
    Min_PSDR = zeros(1, num_gms);                 % Minimum PSDR across all floors for each GM
    Mean_PSDR = zeros(1, num_gms);                % Mean PSDR across all floors for each GM
    Median_PSDR = zeros(1, num_gms);              % Median PSDR across all floors for each GM
    Percentile75_PSDR = zeros(1, num_gms);        % 75th percentile PSDR across all floors for each GM
    Percentile25_PSDR = zeros(1, num_gms);        % 25th percentile PSDR across all floors for each GM
    Percentile65_PSDR = zeros(1, num_gms);        % 65th percentile PSDR across all floors for each GM
    Percentile35_PSDR = zeros(1, num_gms);        % 35th percentile PSDR across all floors for each GM
    
    % Calculate PSDR statistics for each ground motion
    for gm = 1:num_gms
        psdr_values = abs(PSDR(:, gm));           % Get absolute PSDR values for all floors for current GM
        Max_PSDR(1, gm) = max(psdr_values);      % Maximum PSDR across all floors
        Mean_PSDR(1, gm) = mean(psdr_values);    % Average PSDR across all floors
        Median_PSDR(1, gm) = median(psdr_values); % Median PSDR across all floors
        Percentile75_PSDR(1, gm) = prctile(psdr_values, 75); % 75th percentile (upper quartile)
        Percentile25_PSDR(1, gm) = prctile(psdr_values, 25); % 25th percentile (lower quartile)
        Percentile65_PSDR(1, gm) = prctile(psdr_values, 65); % 65th percentile
        Percentile35_PSDR(1, gm) = prctile(psdr_values, 35); % 35th percentile
        Min_PSDR(gm) = min(psdr_values);          % Minimum PSDR across all floors
    end
    
    % Initialize statistical arrays for Peak Floor Acceleration (PFA)
    Max_PFA = zeros(1, num_gms);                  % Maximum PFA across all floors for each GM
    Min_PFA = zeros(1, num_gms);                  % Minimum PFA across all floors for each GM
    Mean_PFA = zeros(1, num_gms);                 % Mean PFA across all floors for each GM
    Median_PFA = zeros(1, num_gms);               % Median PFA across all floors for each GM
    Percentile75_PFA = zeros(1, num_gms);         % 75th percentile PFA across all floors for each GM
    Percentile25_PFA = zeros(1, num_gms);         % 25th percentile PFA across all floors for each GM
    Percentile65_PFA = zeros(1, num_gms);         % 65th percentile PFA across all floors for each GM
    Percentile35_PFA = zeros(1, num_gms);         % 35th percentile PFA across all floors for each GM
    
    % Calculate PFA statistics for each ground motion
    for gm = 1:num_gms
        pfa_values = abs(PFA(:, gm));             % Get absolute PFA values for all floors for current GM
        Max_PFA(1, gm) = max(pfa_values);        % Maximum PFA across all floors
        Mean_PFA(1, gm) = mean(pfa_values);      % Average PFA across all floors
        Median_PFA(1, gm) = median(pfa_values);  % Median PFA across all floors
        Percentile75_PFA(1, gm) = prctile(pfa_values, 75); % 75th percentile (upper quartile)
        Percentile25_PFA(1, gm) = prctile(pfa_values, 25); % 25th percentile (lower quartile)
        Percentile65_PFA(1, gm) = prctile(pfa_values, 65); % 65th percentile
        Percentile35_PFA(1, gm) = prctile(pfa_values, 35); % 35th percentile
        Min_PFA(gm) = min(pfa_values);            % Minimum PFA across all floors
    end
    
    % Initialize statistical arrays for Residual Story Drift Ratio (RSDR)
    Max_RSDR = zeros(1, num_gms);                 % Maximum RSDR across all floors for each GM
    Mean_RSDR = zeros(1, num_gms);                % Mean RSDR across all floors for each GM
    Median_RSDR = zeros(1, num_gms);              % Median RSDR across all floors for each GM
    
    % Calculate RSDR statistics for each ground motion
    for gm = 1:num_gms
        rsdr_values = abs(RSDR(:, gm));           % Get absolute RSDR values for all floors for current GM
        Max_RSDR(1, gm) = max(rsdr_values);      % Maximum RSDR across all floors
        Mean_RSDR(1, gm) = mean(rsdr_values);    % Average RSDR across all floors
        Median_RSDR(1, gm) = median(rsdr_values); % Median RSDR across all floors
    end
    
    %% Define Integration Ranges and Load Fragility Data
    % ---------------------------------------------------------------------
    % Set up integration ranges for loss calculations and load component fragility/cost data
    % ---------------------------------------------------------------------
    
    % Integration ranges for loss calculation (numerical integration limits)
    x_PIDR = (0.0002:0.0002:3)';                  % PIDR range for integration [rad] - fine discretization for accuracy
    x_PFA = (0.02:0.02:10)';                      % PFA range for integration [g] - appropriate range for accelerations
    
    % Load fragility and cost data from HAZUS methodology files
    FragilityDataStructural = readmatrix('HAZUS\Fragility_Data_Structural.csv');
    FragilityDataNonstructuralDrift = readmatrix('HAZUS\Fragility_Data_NonstructuralDrift.csv');
    FragilityDataNonstructuralAcc = readmatrix('HAZUS\Fragility_Data_NonstructuralAcc.csv');
    CostData = readmatrix('HAZUS\Cost_Data.csv');
    
    % Extract median demand thresholds for different damage states
    % Columns contain: [DS1_median, DS1_beta, DS2_median, DS2_beta, DS3_median, DS3_beta, DS4_median, DS4_beta]
    mediandemandStructural = [FragilityDataStructural(:,1), FragilityDataStructural(:,3), ...
                             FragilityDataStructural(:,5), FragilityDataStructural(:,7)];
    mediandemandNonstructuralDrift = [FragilityDataNonstructuralDrift(:,1), FragilityDataNonstructuralDrift(:,3), ...
                                     FragilityDataNonstructuralDrift(:,5), FragilityDataNonstructuralDrift(:,7)];
    mediandemandNonstructuralAcc = [FragilityDataNonstructuralAcc(:,1), FragilityDataNonstructuralAcc(:,3), ...
                                   FragilityDataNonstructuralAcc(:,5), FragilityDataNonstructuralAcc(:,7)];
    
    % Dispersion parameters for fragility functions (logarithmic standard deviations)
    disp_Structural = 0.4;                       % Structural component fragility dispersion
    disp_NonstructuralDrift = 0.5;               % Non-structural drift-sensitive component dispersion
    disp_NonstructuralAcc = 0.6;                 % Non-structural acceleration-sensitive component dispersion
    
    % Calculate coefficient of variation (not used in current analysis but available)
    COV_Structural = disp_Structural ./ mediandemandStructural;
    COV_NonstructuralDrift = disp_NonstructuralDrift ./ mediandemandNonstructuralDrift;
    COV_NonstructuralAcc = disp_NonstructuralAcc ./ mediandemandNonstructuralAcc;
    
    %% Building-Specific Parameters
    % ---------------------------------------------------------------------
    % Extract building geometry and set component-specific parameters based on building height
    % ---------------------------------------------------------------------
    
    % Extract building geometry from loaded data
    Num_stories = Geometry(1, 1);                % Number of stories in building
    Bay_width = Geometry(1, 6);                  % Bay width in feet
    
    % Set median EDP thresholds based on building height category
    % Different building heights have different fragility parameters
    if Num_stories == 1
        % Low-rise building (1 story) - use row 1 of fragility data
        medianEDP_Struct = mediandemandStructural(1, :) / 216;           % Convert to IDR by dividing by conversion factor
        medianEDP_NSDrift = mediandemandNonstructuralDrift(1, :) / 216;  % Convert to IDR
        medianEDP_NSAcc = mediandemandNonstructuralAcc(1, :);            % PFA values used directly (no conversion)
    elseif Num_stories == 5
        % Mid-rise building (5 stories) - use row 2 of fragility data
        medianEDP_Struct = mediandemandStructural(2, :) / 540;           % Convert to IDR with different factor
        medianEDP_NSDrift = mediandemandNonstructuralDrift(2, :) / 540;  % Convert to IDR
        medianEDP_NSAcc = mediandemandNonstructuralAcc(2, :);            % PFA values used directly
    elseif Num_stories > 5
        % High-rise building (>5 stories) - use row 3 of fragility data
        medianEDP_Struct = mediandemandStructural(3, :) / 1123;          % Convert to IDR with high-rise factor
        medianEDP_NSDrift = mediandemandNonstructuralDrift(3, :) / 1123; % Convert to IDR
        medianEDP_NSAcc = mediandemandNonstructuralAcc(3, :);            % PFA values used directly
    end
    
    % Component cost data (repair cost ratios for different damage states)
    CostData_Struct = CostData(1, :);            % Structural component costs (row 1)
    CostData_NSDrift = CostData(3, :);           % Non-structural drift-sensitive costs (row 3)
    CostData_NSAcc = CostData(2, :);             % Non-structural acceleration-sensitive costs (row 2)
    
    %% Alpha Calculation Using Regression Model
    % ---------------------------------------------------------------------
    % Determine alpha value for current building based on fundamental period
    % Uses pre-calibrated regression equation from optimization analysis
    % ---------------------------------------------------------------------
    
    % Extract fundamental period for current building from building information file
    TimePeriod = csvFilePath1(bd+1, 19);         % Period is in column 19, row bd+1 (accounting for header)
    
    % Regression equations for different EDP statistic combinations
    % (Commented equations show alternative regression models based on regression between alpha values and time period)
    
    % Linear fit for MAX and MEAN combination:
    % log_alpha = -1.4475*(log(TimePeriod))+0.751;
    
    % Polynomial (degree 3) fit for MAX and MEAN combination:
    % log_alpha = 1.168*(log(TimePeriod)^3)-0.5263*(log(TimePeriod)^2)-1.9611*(log(TimePeriod))+0.9254;
    
    % Current regression model: Linear fit for 75th and 25th percentile combination:
    log_alpha = -0.20452*(log(TimePeriod))+1.0144;
    
    % Linear fit for MIN and MAX combination:
    % log_alpha = -0.9743*(log(TimePeriod))+0.5204;
    
    % % Current regression model: Linear fit for 35th and 65th percentile combination
    % log_alpha = -0.20187*(log(TimePeriod))+1.4919;
    
    % Convert from log scale to actual alpha value
    alpha = exp(log_alpha);                      % Alpha value for current building
    
    fprintf('  Period: %.3f s, Alpha: %.4f\n', TimePeriod, alpha);
    
    %% Calculate Weighted EDP Statistics
    % ---------------------------------------------------------------------
    % Compute weighted averages using regression-determined alpha value
    % Different weighting schemes available (most are commented out)
    % ---------------------------------------------------------------------
    
    % Alternative weighting schemes (commented out):
    
    % Max-Mean combination:
    % Weighted_PSDR = alpha * Max_PSDR + (1 - alpha) * Mean_PSDR;    
    % Weighted_PFA = alpha * Max_PFA + (1 - alpha) * Mean_PFA;
    
    % 75th-25th percentile combination:
    Weighted_PSDR = alpha * Percentile75_PSDR + (1 - alpha) * Percentile25_PSDR;
    Weighted_PFA = alpha * Percentile75_PFA + (1 - alpha) * Percentile25_PFA;
    
    % Min-Max combination:
    % Weighted_PSDR = alpha * Max_PSDR + (1 - alpha) * Min_PSDR;
    % Weighted_PFA = alpha * Max_PFA + (1 - alpha) * Min_PFA;
    
    % Current approach: 65th-35th percentile combination
    % This provides a balance between conservative (65th) and moderate (35th) estimates
    % Weighted_PSDR = alpha * Percentile65_PSDR + (1 - alpha) * Percentile35_PSDR;
    % Weighted_PFA = alpha * Percentile65_PFA + (1 - alpha) * Percentile35_PFA;
    
    %% Structural Component Loss Assessment
    % ---------------------------------------------------------------------
    % Calculate structural losses using weighted PSDR
    % This follows the HAZUS assembly-based methodology
    % ---------------------------------------------------------------------
    
    % Develop demand model for weighted PSDR using log-linear regression
    % log(EDP) = a + b*log(IM) + error, where IM is intensity measure (Sa)
    md_PIDR = fitlm(log(IM_GM(bd, :)), log(Weighted_PSDR(1, :)));
    PIDR_DemandPar(1:2) = md_PIDR.Coefficients{1:2, 1};  % Extract [intercept, slope]
    PIDR_DemandPar(3) = md_PIDR.RMSE;                     % Extract model uncertainty (sigma)
    SIGMA = PIDR_DemandPar(3);                            % Standard deviation of log-residuals
    
    % Initialize loss arrays for structural components
    NormLoss_Struc = zeros(1, length(SaCalc));           % Normalized structural loss for each Sa level
    NormLoss_Struc_edpj = zeros(length(x_PIDR), 1);      % Pre-allocate temporary array for EDP integration
    
    % Calculate structural losses for each spectral acceleration level
    for im = 1:length(SaCalc)
        % Calculate mean of log(PIDR) given current Sa level using demand model
        MD_PIDR = PIDR_DemandPar(2) * log(SaCalc(im)) + PIDR_DemandPar(1);
        
        % Calculate probability density function of PIDR given Sa
        PDF_PIDR_IM = lognpdf(x_PIDR, MD_PIDR, SIGMA);
        
        % Initialize fragility and loss arrays for current Sa level
        Frag_DS = zeros(length(x_PIDR), length(medianEDP_Struct));     % Fragility for each damage state
        Prob_ds = zeros(length(x_PIDR), length(medianEDP_Struct));     % Probability of being in each damage state
        Loss_EDP = zeros(length(x_PIDR), length(medianEDP_Struct));    % Loss for each damage state
        
        % Calculate fragility functions for each damage state
        for ds = 1:length(medianEDP_Struct)
            % Lognormal CDF: P(DS ≥ ds | EDP) = Φ[(ln(EDP) - ln(median))/β]
            Frag_DS(:, ds) = logncdf(x_PIDR, log(medianEDP_Struct(1, ds)), disp_Structural);
        end
        
        % Calculate probability of being in each discrete damage state
        % P(DS = ds) = P(DS ≥ ds) - P(DS ≥ ds+1)
        for p = 1:length(medianEDP_Struct)-1
            Prob_ds(:, p) = Frag_DS(:, p) - Frag_DS(:, p+1);  % Probability of damage state p
            Prob_ds(:, p+1) = Frag_DS(:, p+1);                % Probability of highest damage state
        end
        
        % Calculate expected loss for each damage state
        for p = 1:length(medianEDP_Struct)
            % Expected loss = Cost ratio × Probability of DS × PDF of EDP given IM
            Loss_EDP(:, p) = CostData_Struct(1, p) * (Prob_ds(:, p) .* PDF_PIDR_IM);
        end
        
        % Sum across all damage states and integrate over EDP range
        NormLoss_Struc_edpj = sum(Loss_EDP, 2);              % Sum over damage states
        NormLoss_Struc(1, im) = trapz(x_PIDR, NormLoss_Struc_edpj); % Integrate over PIDR range
    end
    
    %% Non-Structural Drift-Sensitive Loss Assessment
    % ---------------------------------------------------------------------
    % Calculate drift-sensitive non-structural losses using weighted PSDR
    % Similar methodology to structural components but different fragility/cost data
    % ---------------------------------------------------------------------
    
    % Initialize loss arrays for non-structural drift-sensitive components
    NormLoss_NonStrucDrift = zeros(1, length(SaCalc));      % Normalized loss for each Sa level
    NormLoss_NonStrucDrift_edpj = zeros(length(x_PIDR), 1); % Pre-allocate temporary array
    
    % Calculate non-structural drift-sensitive losses for each Sa level
    for im = 1:length(SaCalc)
        % Use same demand model as structural (both driven by drift)
        MD_PIDR = PIDR_DemandPar(2) * log(SaCalc(im)) + PIDR_DemandPar(1);
        PDF_PIDR_IM = lognpdf(x_PIDR, MD_PIDR, SIGMA);
        
        % Initialize arrays for non-structural drift-sensitive components
        Frag_DS_NSD = zeros(length(x_PIDR), length(medianEDP_NSDrift));
        Prob_ds_NSD = zeros(length(x_PIDR), length(medianEDP_NSDrift));
        Loss_EDP_NSD = zeros(length(x_PIDR), length(medianEDP_NSDrift));
        
        % Calculate fragility functions for non-structural drift-sensitive components
        for ds = 1:length(medianEDP_NSDrift)
            Frag_DS_NSD(:, ds) = logncdf(x_PIDR, log(medianEDP_NSDrift(1, ds)), disp_NonstructuralDrift);
        end
        
        % Calculate probability of being in each damage state
        for p = 1:length(medianEDP_NSDrift)-1
            Prob_ds_NSD(:, p) = Frag_DS_NSD(:, p) - Frag_DS_NSD(:, p+1);
            Prob_ds_NSD(:, p+1) = Frag_DS_NSD(:, p+1);
        end
        
        % Calculate expected loss for each damage state
        for p = 1:length(medianEDP_NSDrift)
            Loss_EDP_NSD(:, p) = CostData_NSDrift(1, p) * (Prob_ds_NSD(:, p) .* PDF_PIDR_IM);
        end
        
        % Sum and integrate to get total loss
        NormLoss_NonStrucDrift_edpj = sum(Loss_EDP_NSD, 2);
        NormLoss_NonStrucDrift(1, im) = trapz(x_PIDR, NormLoss_NonStrucDrift_edpj);
    end
    
    %% Non-Structural Acceleration-Sensitive Loss Assessment
    % ---------------------------------------------------------------------
    % Calculate acceleration-sensitive non-structural losses using weighted PFA
    % Uses acceleration demand model instead of drift
    % ---------------------------------------------------------------------
    
    % Develop demand model for weighted PFA using log-linear regression
    md_PFA = fitlm(log(IM_GM(bd, :)), log(Weighted_PFA(1, :)));
    PFA_DemandPar(1:2) = md_PFA.Coefficients{1:2, 1};    % Extract regression coefficients
    PFA_DemandPar(3) = md_PFA.RMSE;                       % Extract model uncertainty
    SIGMA_PFA = PFA_DemandPar(3);                         % Standard deviation for PFA model
    
    % Initialize loss arrays for non-structural acceleration-sensitive components
    NormLoss_NonStrucAcc = zeros(1, length(SaCalc));      % Normalized loss for each Sa level
    NormLoss_NonStrucAcc_edpj = zeros(length(x_PFA), 1);  % Pre-allocate temporary array
    
    % Calculate non-structural acceleration-sensitive losses for each Sa level
    for im = 1:length(SaCalc)
        % Calculate mean of log(PFA) given current Sa level
        MD_PFA = PFA_DemandPar(2) * log(SaCalc(im)) + PFA_DemandPar(1);
        
        % Calculate probability density function of PFA given Sa
        PDF_PFA_IM = lognpdf(x_PFA, MD_PFA, SIGMA_PFA);
        
        % Initialize arrays for non-structural acceleration-sensitive components
        Frag_DS_NSA = zeros(length(x_PFA), length(medianEDP_NSAcc));
        Prob_ds_NSA = zeros(length(x_PFA), length(medianEDP_NSAcc));
        Loss_EDP_NSA = zeros(length(x_PFA), length(medianEDP_NSAcc));
        
        % Calculate fragility functions for acceleration-sensitive components
        for ds = 1:length(medianEDP_NSAcc)
            Frag_DS_NSA(:, ds) = logncdf(x_PFA, log(medianEDP_NSAcc(1, ds)), disp_NonstructuralAcc);
        end
        
        % Calculate probability of being in each damage state
        for p = 1:length(medianEDP_NSAcc)-1
            Prob_ds_NSA(:, p) = Frag_DS_NSA(:, p) - Frag_DS_NSA(:, p+1);
            Prob_ds_NSA(:, p+1) = Frag_DS_NSA(:, p+1);
        end
        
        % Calculate expected loss for each damage state
        for p = 1:length(medianEDP_NSAcc)
            Loss_EDP_NSA(:, p) = CostData_NSAcc(1, p) * (Prob_ds_NSA(:, p) .* PDF_PFA_IM);
        end
        
        % Sum and integrate to get total loss
        NormLoss_NonStrucAcc_edpj = sum(Loss_EDP_NSA, 2);
        NormLoss_NonStrucAcc(1, im) = trapz(x_PFA, NormLoss_NonStrucAcc_edpj);
    end
    
    %% Combine Component Losses
    % ---------------------------------------------------------------------
    % Combine all component losses (excluding demolition effects)
    % Convert from percentage to fraction and sum all components
    % ---------------------------------------------------------------------
    
    % Convert normalized losses from percentage to fraction of replacement cost
    NormLoss_Struc = NormLoss_Struc(1, :) / 100;           % Convert from % to fraction
    NormLoss_NonStrucDrift = NormLoss_NonStrucDrift(1, :) / 100; % Convert from % to fraction
    NormLoss_NonStrucAcc = NormLoss_NonStrucAcc(1, :) / 100;     % Convert from % to fraction
    
    % Calculate total repair loss (before considering demolition)
    NormTotLoss = NormLoss_Struc + NormLoss_NonStrucDrift + NormLoss_NonStrucAcc;
    
    %% Demolition Loss Assessment
    % ---------------------------------------------------------------------
    % Calculate demolition probability and associated losses
    % Based on residual drift exceeding acceptable thresholds (Ramirez & Miranda 2012)
    % ---------------------------------------------------------------------
    
    % Define integration range for residual story drift ratio
    x_RSDR_pdf = (0.0002:0.0002:0.06)';                    % RSDR range for integration [rad]
    
    % Calculate demolition fragility function based on Ramirez and Miranda (2012)
    % P(Demolition | RSDR) = Φ[(ln(RSDR) - ln(θ))/β]
    PD_RSDR = normcdf((log(x_RSDR_pdf / theta_hat_Demolish)) / beta_hat_Demolish);
    
    % Initialize demolition probability array for each Sa level
    PD_IM_NC = zeros(length(SaCalc), 1);                   % Probability of demolition given Sa
    
    % Calculate demolition probability for each spectral acceleration level
    for im = 1:length(SaCalc)
        
        % Develop RSDR demand model using maximum RSDR (most critical for demolition decisions)
        md_RSDR = fitlm(log(IM_GM(bd, :)), log(Max_RSDR(1, :)));
        
        % Extract regression coefficients and model uncertainty
        RSDR_DemandPar(1:2) = md_RSDR.Coefficients{1:2, 1};   % [intercept, slope]
        RSDR_DemandPar(3) = md_RSDR.RMSE;                      % Model uncertainty (sigma)
        SIGMA_RSDR = RSDR_DemandPar(3);                        % Standard deviation for RSDR model
        
        % Calculate mean of log(RSDR) given current Sa level
        MD_RSDR = RSDR_DemandPar(2) * log(SaCalc(im)) + RSDR_DemandPar(1);
        
        % Calculate probability density function of RSDR given Sa
        dP_RIDR_IM = lognpdf(x_RSDR_pdf, MD_RSDR, SIGMA_RSDR);
        
        % Convolve demolition fragility with RSDR distribution
        PD_RSDR_x_dP_RSDR_IM = PD_RSDR(:, 1) .* dP_RIDR_IM(:, 1);
        
        % Integrate to get total demolition probability at this Sa level
        PD_IM_NC(im, 1) = trapz(x_RSDR_pdf, PD_RSDR_x_dP_RSDR_IM);
    end
    
    % Transpose for consistency with other arrays
    PD_IM_NC = PD_IM_NC';
    
    %% Calculate Total Costs and Final Losses
    % ---------------------------------------------------------------------
    % Compute replacement costs and apply demolition effects to all loss components
    % ---------------------------------------------------------------------
    
    % Calculate building replacement cost based on geometry
    Num_bays = 5;                                          % Assumed number of bays in each direction
    Floor_area = (Bay_width * Num_bays) * (Bay_width * Num_bays); % Total floor area [ft²]
    Per_sqft_Cost = 250;                                   % Construction cost per square foot [$/ft²]
    Tot_Replac_Cost = Num_stories * Floor_area * Per_sqft_Cost;   % Total replacement cost [$]
    
    % Calculate demolition-related costs
    DemoLoss = 0.1 * Tot_Replac_Cost;                      % Demolition cost (10% of replacement cost)
    Tot_Loss_Demo = Tot_Replac_Cost + DemoLoss;            % Total loss if building is demolished
    
    % Apply demolition probability to repair losses
    % If building is demolished, no repair is performed, so multiply repair losses by (1 - P_demolition)
    Norm_Loss_Struc(1, :) = NormLoss_Struc(1, :) .* (1 - PD_IM_NC(1, :));
    Loss_Struc(1, :) = Norm_Loss_Struc(1, :) .* Tot_Replac_Cost;
    
    Norm_Loss_NonStruc_Drift(1, :) = NormLoss_NonStrucDrift(1, :) .* (1 - PD_IM_NC(1, :));
    Loss_NonStruc_Drift(1, :) = Norm_Loss_NonStruc_Drift(1, :) .* Tot_Replac_Cost;
    
    Norm_Loss_NonStruc_Acc(1, :) = NormLoss_NonStrucAcc(1, :) .* (1 - PD_IM_NC(1, :));
    Loss_NonStruc_Acc(1, :) = Norm_Loss_NonStruc_Acc(1, :) .* Tot_Replac_Cost;
    
    % Calculate demolition loss (occurs with probability P_demolition)
    Loss_Demo(1, :) = Tot_Loss_Demo * PD_IM_NC(1, :);
    Norm_Loss_Demo(1, :) = Loss_Demo(1, :) ./ Tot_Replac_Cost;  % Normalize by replacement cost
    
    % Calculate total normalized loss (repair OR demolition, mutually exclusive)
    Norm_Tot_Loss(1, :) = NormTotLoss(1, :) .* (1 - PD_IM_NC(1, :)) + Norm_Loss_Demo;
    Loss_Tot(1, :) = Norm_Tot_Loss(1, :) .* Tot_Replac_Cost;
    
    %% Seismic Hazard Processing and EAL Calculation
    % ---------------------------------------------------------------------
    % Load hazard data and calculate Expected Annual Loss through integration
    % ---------------------------------------------------------------------
    
    % Load seismic hazard data from USGS
    Hazard_file = 'USGSHazard_data.csv';
    HazardData = readmatrix(Hazard_file);
    AllLamda_IM = HazardData(:, 2:12);                     % Hazard rates for different periods [1/year]
    IM_values = HazardData(:, 1);                          % Sa values for hazard curve [g]
    
    % Extract hazard rates for different structural periods
    a = AllLamda_IM(:, 2);                                 % T = 0.1s hazard rates
    b = AllLamda_IM(:, 3);                                 % T = 0.2s hazard rates
    c = AllLamda_IM(:, 4);                                 % T = 0.3s hazard rates
    d = AllLamda_IM(:, 5);                                 % T = 0.5s hazard rates
    e = AllLamda_IM(:, 6);                                 % T = 0.75s hazard rates
    f = AllLamda_IM(:, 7);                                 % T = 1.0s hazard rates
    g = AllLamda_IM(:, 8);                                 % T = 2.0s hazard rates
    h = AllLamda_IM(:, 9);                                 % T = 3.0s hazard rates
    i = AllLamda_IM(:, 10);                                % T = 4.0s hazard rates
    j = AllLamda_IM(:, 11);                                % T = 5.0s hazard rates
    
    % Interpolate hazard curve for building's specific fundamental period
    T_target = Target_time_period(bd, 1);                  % Current building's fundamental period
    Lamda_IM_target = zeros(length(IM_values), length(Target_time_period)); % Pre-allocate
    
    % Period-based interpolation of hazard curves (logarithmic interpolation)
    if T_target >= 0.1 && T_target < 0.2
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([0.1, 0.2], [log(a(nn, 1)), log(b(nn, 1))], T_target, 'linear');
        end
    elseif T_target >= 0.2 && T_target < 0.3
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([0.2, 0.3], [log(b(nn, 1)), log(c(nn, 1))], T_target, 'linear');
        end
    elseif T_target >= 0.3 && T_target < 0.5
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([0.3, 0.5], [log(c(nn, 1)), log(d(nn, 1))], T_target, 'linear');
        end
    elseif T_target >= 0.5 && T_target < 0.75
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([0.5, 0.75], [log(d(nn, 1)), log(e(nn, 1))], T_target, 'linear');
        end
    elseif T_target >= 0.75 && T_target < 1
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([0.75, 1], [log(e(nn, 1)), log(f(nn, 1))], T_target, 'linear');
        end
    elseif T_target >= 1 && T_target < 2
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([1, 2], [log(f(nn, 1)), log(g(nn, 1))], T_target, 'linear');
        end
    elseif T_target >= 2 && T_target < 3
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([2, 3], [log(g(nn, 1)), log(h(nn, 1))], T_target, 'linear');
        end
    elseif T_target >= 3 && T_target < 4
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([3, 4], [log(h(nn, 1)), log(i(nn, 1))], T_target, 'linear');
        end
    elseif T_target >= 4 && T_target < 5
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([4, 5], [log(i(nn, 1)), log(j(nn, 1))], T_target, 'linear');
        end
    end
    
    % Convert interpolated hazard curve to exceedance rates
    hazardCurveSa = [IM_values, exp(Lamda_IM_target(:, bd))];
    
    % Numerical differentiation of hazard curve to get hazard density function
    % This converts from exceedance rate λ(Sa) to probability density d λ/d Sa
    PDF_temp = zeros(length(SaCalc), 1);                   % Pre-allocate hazard density array
    hh = 0.0000001;                                        % Small increment for finite difference
    
    for kk = 1:length(SaCalc)
        % Two-point finite difference approximation: f'(x) ≈ [f(x+h) - f(x-h)]/(2h)
        fxp = interp1(hazardCurveSa(:, 1), hazardCurveSa(:, 2), SaCalc(kk) + hh, 'spline');
        fxn = interp1(hazardCurveSa(:, 1), hazardCurveSa(:, 2), SaCalc(kk) - hh, 'spline');
        PDF_temp(kk) = abs((fxp - fxn)) / (2 * hh);        % Absolute value of derivative
    end
    
    % Smooth the hazard density function to reduce numerical noise
    dlambdaSa(:, bd) = smooth(PDF_temp);
    
    %% Expected Annual Loss (EAL) Computation
    % ---------------------------------------------------------------------
    % Calculate EAL by integrating loss × hazard density for all loss components
    % EAL = ∫ Loss(Sa) × |dλ/dSa| dSa
    % ---------------------------------------------------------------------
    
    % Convolution of losses with hazard density function
    Conv_loss_Struc = Loss_Struc(1, :)' .* dlambdaSa(:, bd);           % Structural loss convolution
    Conv_loss_NonStrucDrift = Loss_NonStruc_Drift(1, :)' .* dlambdaSa(:, bd); % Non-structural drift loss convolution
    Conv_loss_NonStrucAcc = Loss_NonStruc_Acc(1, :)' .* dlambdaSa(:, bd);    % Non-structural acceleration loss convolution
    Conv_loss_Tot = Loss_Tot(1, :)' .* dlambdaSa(:, bd);               % Total loss convolution
    Conv_loss_Demo = Loss_Demo(1, :)' .* dlambdaSa(:, bd);             % Demolition loss convolution
    
    % Convolution of normalized losses with hazard density function
    Conv_Norm_loss_Struc = Norm_Loss_Struc(1, :)' .* dlambdaSa(:, bd);
    Conv_Norm_loss_NonStrucDrift = Norm_Loss_NonStruc_Drift(1, :)' .* dlambdaSa(:, bd);
    Conv_Norm_loss_NonStrucAcc = Norm_Loss_NonStruc_Acc(1, :)' .* dlambdaSa(:, bd);
    Conv_Norm_loss_Tot = Norm_Tot_Loss(1, :)' .* dlambdaSa(:, bd);
    Conv_Norm_loss_Demo = Norm_Loss_Demo(1, :)' .* dlambdaSa(:, bd);
    
    % Numerical integration using trapezoidal rule to get Expected Annual Loss
    EAL_Tot(bd, 1) = trapz(SaCalc, Conv_loss_Tot);                    % Total EAL [$]
    EAL_Demo(bd, 1) = trapz(SaCalc, Conv_loss_Demo);                  % Demolition EAL [$]
    EAL_Struc(bd, 1) = trapz(SaCalc, Conv_loss_Struc);               % Structural EAL [$]
    EAL_NonStrucDrift(bd, 1) = trapz(SaCalc, Conv_loss_NonStrucDrift); % Non-structural drift EAL [$]
    EAL_NonStrucAcc(bd, 1) = trapz(SaCalc, Conv_loss_NonStrucAcc);    % Non-structural acceleration EAL [$]
    
    % Normalized Expected Annual Loss (with respect to replacement cost)
    Norm_EAL_Tot(bd, 1) = trapz(SaCalc, Conv_Norm_loss_Tot);          % Total normalized EAL [-]
    Norm_EAL_Demo(bd, 1) = trapz(SaCalc, Conv_Norm_loss_Demo);        % Demolition normalized EAL [-]
    Norm_EAL_Struc(bd, 1) = trapz(SaCalc, Conv_Norm_loss_Struc);      % Structural normalized EAL [-]
    Norm_EAL_NonStrucDrift(bd, 1) = trapz(SaCalc, Conv_Norm_loss_NonStrucDrift); % Non-structural drift normalized EAL [-]
    Norm_EAL_NonStrucAcc(bd, 1) = trapz(SaCalc, Conv_Norm_loss_NonStrucAcc);     % Non-structural acceleration normalized EAL [-]
    
    % Display results for current building
    fprintf('  EAL Results: Total=%.0f$, Normalized=%.4f\n', EAL_Tot(bd, 1), Norm_EAL_Tot(bd, 1));
    fprintf('  Building %d completed.\n\n', ID(bd));
    
end  % End main building loop

%% Results Compilation and Output
% =========================================================================
% Compile final results and save to Excel file
% =========================================================================

fprintf('=================================================================\n');
fprintf('LOSS ASSESSMENT COMPLETED SUCCESSFULLY\n');
fprintf('=================================================================\n\n');

% Compile final output matrix with all EAL results
All_Output(:, 1) = ID(:, 1);                              % Building ID numbers
All_Output(:, 2) = EAL_Struc(:, 1);                       % Structural EAL [$]
All_Output(:, 3) = EAL_NonStrucDrift(:, 1);               % Non-structural drift-sensitive EAL [$]
All_Output(:, 4) = EAL_NonStrucAcc(:, 1);                 % Non-structural acceleration-sensitive EAL [$]
All_Output(:, 5) = EAL_Demo(:, 1);                        % Demolition EAL [$]
All_Output(:, 6) = EAL_Tot(:, 1);                         % Total EAL [$]
All_Output(:, 7) = Norm_EAL_Struc(:, 1);                  % Normalized structural EAL [-]
All_Output(:, 8) = Norm_EAL_NonStrucDrift(:, 1);          % Normalized non-structural drift-sensitive EAL [-]
All_Output(:, 9) = Norm_EAL_NonStrucAcc(:, 1);            % Normalized non-structural acceleration-sensitive EAL [-]
All_Output(:, 10) = Norm_EAL_Demo(:, 1);                  % Normalized demolition EAL [-]
All_Output(:, 11) = Norm_EAL_Tot(:, 1);                   % Normalized total EAL [-]

% % Save final results to Excel file (Uncomment below
% lines to save the file)

% excelFileName = 'HAZUS_NormEAL_Weighted_Testing_25_75.xlsx';
% sheetName = 'Sheet1';
% writematrix(All_Output, excelFileName, 'Sheet', sheetName);
% 
% fprintf('Results saved to: %s\n', excelFileName);


fprintf('\n=================================================================\n');
fprintf('ANALYSIS COMPLETE - RESULTS READY FOR POST-PROCESSING\n');
fprintf('=================================================================\n');

%% End of Script