%% HAZUS Seismic Loss Assessment Main Code 
% =========================================================================
% HAZUS-BASED SEISMIC LOSS ASSESSMENT FOR STEEL MOMENT-RESISTING FRAMES
% =========================================================================
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
% Description:
%   This script implements the HAZUS methodology for seismic loss assessment
%   of steel moment-resisting frame buildings. HAZUS uses simplified fragility
%   and cost functions for different assemblies (performance groups) compared
%   to the more detailed component-based approach of FEMA P-58.
%
% COMPARISON: HAZUS vs FEMA P-58 METHODOLOGIES
% ============================================
% 
% 1. FRAGILITY APPROACH:
%    - HAZUS: Uses building-level fragility curves with 4 damage states
%      (Slight, Moderate, Extensive, Complete) based on structural type
%    - FEMA P-58: Uses component-level fragility curves for individual
%      building components (beams, columns, partitions, equipment, etc.)
%
% 2. LOSS CALCULATION:
%    - HAZUS: Applies normalized loss ratios (% of replacement cost) to
%      each damage state for broad component categories (assemblies or
%      performance groups)
%    - FEMA P-58: Calculates detailed repair costs for each component
%      based on quantities, unit costs, and damage state probabilities
%
% 3. COMPONENT RESOLUTION:
%    - HAZUS: Three main categories: Structural, Non-structural Drift-sensitive,
%      Non-structural Acceleration-sensitive
%    - FEMA P-58: Detailed component inventory (many component types
%      with specific fragility and cost functions)
%
% 4. DEMAND PARAMETERS:
%    - HAZUS: Uses maximum inter-story drift and peak floor acceleration
%      across the entire building
%    - FEMA P-58: Uses floor-by-floor response parameters for detailed
%      spatial distribution of demands
%
% 5. ACCURACY vs EFFICIENCY:
%    - HAZUS: Faster, suitable for regional/portfolio assessment, less precise
%    - FEMA P-58: More accurate, detailed, suitable for individual building
%      assessment, computationally intensive
%
% Reference:
%   FEMA (2003). Multi-hazard Loss Estimation Methodology - Earthquake Model
%   HAZUS®MH MR4 Technical Manual. Federal Emergency Management Agency.
%
% Input Files Required:
%   - GuanDataBase-IMs.csv: Intensity measures for 240 ground motions
%   - Building response files (PSDR, PFA, RSDR) in Guan_database/
%   - HAZUS fragility and cost data files in HAZUS/ directory
%   - USGSHazard_data.csv: Seismic hazard curves
%
% Output:
%   - Expected Annual Loss (EAL) values for each building [$]
%   - Normalized EAL values (as fraction of replacement cost) [-]
%   - Assembly-level breakdown following HAZUS categories

% =========================================================================
% INITIALIZATION AND SETUP
% =========================================================================
clc; 
clearvars; 
close all;
clear;

%% Input Ground Motion Intensity Measures
% Read spectral acceleration values for all buildings and ground motions
IM_FilePath = 'GuanDataBase-IMs.csv';

% Load intensity measure data
Data = readmatrix(IM_FilePath);
IM_GM = Data(:, 3:242);              % IM values for all 240 ground motions (building-wise)
ID = Data(:, 1);                     % Building IDs
IM_GM_trans = IM_GM';                % Transpose for easier processing
Target_time_period = Data(:, 2);     % Fundamental periods of buildings

% Define spectral acceleration range for loss calculations
% Extended range to 6.0g for comprehensive HAZUS assessment
SaCalc = [0.001, 0.01, 0.05, 0.1:0.01:6]'; % Sa values from 0.001g to 6.0g

%% Initialize Output Arrays for Building-wise Assessment
% Pre-allocate arrays for computational efficiency

% Hazard curve derivatives
dlambdaSa = zeros(length(SaCalc), length(ID));

% Expected Annual Loss (EAL) arrays [$]
EAL_Tot = zeros(length(ID), 1);                    % Total EAL
EAL_Demo = zeros(length(ID), 1);                   % Demolition EAL
EAL_Struc = zeros(length(ID), 1);                  % Structural EAL
EAL_NonStrucDrift = zeros(length(ID), 1);          % Non-structural drift-sensitive EAL
EAL_NonStrucAcc = zeros(length(ID), 1);            % Non-structural acceleration-sensitive EAL

% Normalized EAL arrays (fraction of replacement cost) [-]
Norm_EAL_Tot = zeros(length(ID), 1);               % Total normalized EAL
Norm_EAL_Demo = zeros(length(ID), 1);              % Demolition normalized EAL
Norm_EAL_Struc = zeros(length(ID), 1);             % Structural normalized EAL
Norm_EAL_NonStrucDrift = zeros(length(ID), 1);     % Non-structural drift normalized EAL
Norm_EAL_NonStrucAcc = zeros(length(ID), 1);       % Non-structural acceleration normalized EAL

% =========================================================================
% MAIN BUILDING LOOP - PROCESS EACH BUILDING INDIVIDUALLY
% =========================================================================
for bd = 100:length(ID)-520
    
    %% Load Building-Specific Response Data
    % Construct file paths for structural response data
    filePSDR = strcat(['Guan_database\Structural', ...
        'Responses\EDPsUnder240GMs\PeakStoryDrift'], '\', ...
        'Building_', num2str(ID(bd), '%.1d'), '_PSDR.csv');
    
    filePFA = strcat(['Guan_database\Structural', ...
        'Responses\EDPsUnder240GMs\PeakFloorAcceleration'], '\', ...
        'Building_', num2str(ID(bd), '%.1d'), '_PFA.csv');
    
    fileRSDR = strcat(['Guan_database\Structural', ...
        'Responses\EDPsUnder240GMs\ResidualStoryDrift'], '\', ...
        'Building_', num2str(ID(bd), '%.1d'), '_RSDR.csv');
    
    fileGeometry = strcat('Guan_database\BuildingDesigns', '\', ...
        'Building_', num2str(ID(bd), '%.1d'), '\', 'Geometry.csv');

    % =====================================================================
    % RESPONSE DATA PROCESSING - HAZUS APPROACH
    % =====================================================================
    % HAZUS uses maximum values across the building rather than floor-by-floor
    % analysis as in FEMA P-58, representing a simplified assessment approach
    
    % Load and process Peak Story Drift Ratios (PSDR)
    PSDR = readmatrix(filePSDR);        % Inter-story drift ratios (dimensionless)
    
    % Initialize arrays for maximum, mean, and median response parameters
    Max_PSDR = zeros(1, length(PSDR));
    Mean_PSDR = zeros(1, length(PSDR));
    Median_PSDR = zeros(1, length(PSDR));    
    
    % Calculate statistical measures for each ground motion
    for gm = 1:length(PSDR)
        Max_PSDR(1, gm) = max(abs(PSDR(:, gm)));        % Maximum across all stories
        Mean_PSDR(1, gm) = mean(abs(PSDR(:, gm)));       % Average across all stories
        Median_PSDR(1, gm) = median(abs(PSDR(:, gm)));   % Median across all stories
    end    

    % Load and process Peak Floor Accelerations (PFA)
    PFA = readmatrix(filePFA);
    PFA = PFA * 0.0254 * 0.10197;      % Convert from in/s² to g (gravity units)
    
    Max_PFA = zeros(1, length(PFA));
    Mean_PFA = zeros(1, length(PFA));
    Median_PFA = zeros(1, length(PFA));
    
    for gm = 1:length(PFA)
        Max_PFA(1, gm) = max(abs(PFA(:, gm)));          % Maximum across all floors
        Mean_PFA(1, gm) = mean(abs(PFA(:, gm)));         % Average across all floors
        Median_PFA(1, gm) = median(abs(PFA(:, gm)));     % Median across all floors
    end
    
    % Load and process Residual Story Drift Ratios (RSDR)
    RSDR = readmatrix(fileRSDR);       % Permanent drift after earthquake
    
    Max_RSDR = zeros(1, length(RSDR));
    Mean_RSDR = zeros(1, length(RSDR));
    Median_RSDR = zeros(1, length(RSDR));
    
    for gm = 1:length(RSDR)
        Max_RSDR(1, gm) = max(abs(RSDR(:, gm)));        % Maximum residual drift
        Mean_RSDR(1, gm) = mean(abs(RSDR(:, gm)));       % Average residual drift
        Median_RSDR(1, gm) = median(abs(RSDR(:, gm)));   % Median residual drift
    end

    % Define engineering demand parameter ranges for integration
    x_PIDR = (0.0002:0.0002:3)';       % Peak inter-story drift range (dimensionless)
    x_PFA = (0.02:0.02:10)';           % Peak floor acceleration range [g]
  
    % Load building geometry parameters
    Geometry = readmatrix(fileGeometry);

    % =====================================================================
    % HAZUS FRAGILITY AND COST DATA LOADING
    % =====================================================================
    % HAZUS uses standardized fragility curves and loss ratios rather than
    % the detailed component-specific data used in FEMA P-58
    
    % Initialize normalized loss arrays for HAZUS methodology
    NormLoss_Struc = zeros(1, length(SaCalc));
    NormLoss_Struc_edpj = zeros(length(x_PIDR), 1);
    NormLoss_NonStrucDrift = zeros(1, length(SaCalc));
    NormLoss_NonStrucDrift_edpj = zeros(length(x_PIDR), 1);
    NormLoss_NonStrucAcc = zeros(1, length(SaCalc));
    NormLoss_NonStrucAcc_edpj = zeros(length(x_PFA), 1);

    % Load HAZUS fragility and cost data files
    fileFragility_Data_Structural = strcat('HAZUS', '\', 'Fragility_Data_Structural.csv');
    fileFragility_Data_NonstructuralDrift = strcat('HAZUS', '\', 'Fragility_Data_NonstructuralDrift.csv');    
    fileFragility_Data_NonstructuralAcc = strcat('HAZUS', '\', 'Fragility_Data_NonstructuralAcc.csv');    
    fileCost_Data = strcat('HAZUS', '\', 'Cost_Data.csv');    
    
    % Read HAZUS data matrices
    FragilityDataStructural = readmatrix(fileFragility_Data_Structural);
    FragilityDataNonstructuralDrift = readmatrix(fileFragility_Data_NonstructuralDrift);
    FragilityDataNonstructuralAcc = readmatrix(fileFragility_Data_NonstructuralAcc);
    CostData = readmatrix(fileCost_Data);
    
    % Extract median demand values in terms of spectral displacement
    % HAZUS fragility curves are organized by damage state (Slight, Moderate, Extensive, Complete)
    mediandemandStructural = [FragilityDataStructural(:, 1), FragilityDataStructural(:, 3), ...
                             FragilityDataStructural(:, 5), FragilityDataStructural(:, 7)];
    
    mediandemandNonstructuralDrift = [FragilityDataNonstructuralDrift(:, 1), FragilityDataNonstructuralDrift(:, 3), ...
                                     FragilityDataNonstructuralDrift(:, 5), FragilityDataNonstructuralDrift(:, 7)];
    
    mediandemandNonstructuralAcc = [FragilityDataNonstructuralAcc(:, 1), FragilityDataNonstructuralAcc(:, 3), ...
                                   FragilityDataNonstructuralAcc(:, 5), FragilityDataNonstructuralAcc(:, 7)];

    % HAZUS uses simplified dispersion values rather than component-specific dispersions
    disp_Structural = 0.4;             % Log-standard deviation for structural components
    disp_NonstructuralDrift = 0.5;     % Log-standard deviation for drift-sensitive components
    disp_NonstructuralAcc = 0.6;       % Log-standard deviation for acceleration-sensitive components

    % Calculate coefficients of variation (not used in current implementation)
    COV_Structural = disp_Structural ./ mediandemandStructural;
    COV_NonstructuralDrift = disp_NonstructuralDrift ./ mediandemandNonstructuralDrift;
    COV_NonstructuralAcc = disp_NonstructuralAcc ./ mediandemandNonstructuralAcc;

    % =====================================================================
    % STRUCTURAL COMPONENT LOSS ASSESSMENT - HAZUS APPROACH
    % =====================================================================
    % HAZUS uses building-level fragility curves rather than individual
    % component assessment as in FEMA P-58
    
    % Extract building parameters
    Num_stories = Geometry(1, 1);      % Number of stories
    Bay_width = Geometry(1, 6);        % Bay width [ft]
    
    % Select appropriate fragility parameters based on building height
    % HAZUS categorizes buildings by height ranges with different vulnerabilities
    if Num_stories == 1
        % Low-rise buildings: Convert spectral displacement to inter-story drift
        medianEDP_Struct = mediandemandStructural(1, :) / 216;   % Converted to IDR
        
    elseif Num_stories == 5
        % Mid-rise buildings
        medianEDP_Struct = mediandemandStructural(2, :) / 540;   % Converted to IDR
        
    elseif Num_stories > 5
        % High-rise buildings
        medianEDP_Struct = mediandemandStructural(3, :) / 1123;  % Converted to IDR
    end  
    
    % Use standardized dispersion value for structural components
    dispEDP_Struct = disp_Structural;  
    
    % Extract cost ratios for structural components (% of replacement cost)
    CostData_Struct = CostData(1, :);

    % Develop demand model: Sa vs Maximum Inter-story Drift Ratio
    % This represents the HAZUS approach of using building-level maximum response
    md_PIDR = fitlm(log(IM_GM(bd, :)), log(Max_PSDR(1, :)));
    PIDR_DemandPar(1:2) = md_PIDR.Coefficients{1:2, 1};    % [intercept, slope]
    PIDR_DemandPar(3) = md_PIDR.RMSE;                       % Residual standard error
    SIGMA = PIDR_DemandPar(3);                              % Log-standard deviation

    % Loop through each intensity measure level
    for im = 1:length(SaCalc)

        % Calculate mean log demand for current intensity level
        MD_PIDR = PIDR_DemandPar(2) * log(SaCalc(im)) + PIDR_DemandPar(1);
        
        % Calculate probability density function of drift given Sa
        PDF_PIDR_IM = lognpdf(x_PIDR, MD_PIDR, SIGMA);
        
        % Initialize damage state calculation arrays
        Frag_DS = zeros(length(x_PIDR), length(medianEDP_Struct));
        Prob_ds = zeros(length(x_PIDR), length(medianEDP_Struct));
        Loss_EDP = zeros(length(x_PIDR), length(medianEDP_Struct));                                      
    
        % Calculate fragility functions for each damage state
        for ds = 1:length(medianEDP_Struct)
            % Lognormal CDF for exceedance probability P(DS ≥ ds | EDP)
            Frag_DS(:, ds) = logncdf(x_PIDR, log(medianEDP_Struct(1, ds)), dispEDP_Struct);
        end               

        % Calculate probability of being in each damage state
        % P(DS = i) = P(DS ≥ i) - P(DS ≥ i+1)
        for p = 1:length(medianEDP_Struct)-1
            Prob_ds(:, p) = Frag_DS(:, p) - Frag_DS(:, p+1);
            Prob_ds(:, p+1) = Frag_DS(:, p+1);
        end
    
        % Calculate expected loss for each damage state
        for p = 1:length(medianEDP_Struct)
            % Loss = Cost_Ratio × P(DS) × PDF(EDP|IM)
            Loss_EDP(:, p) = CostData_Struct(1, p) * (Prob_ds(:, p) .* PDF_PIDR_IM); 
        end
        
        % Sum losses across all damage states
        NormLoss_Struc_edpj = sum(Loss_EDP, 2);
        
        % Integrate over demand range to get expected loss
        NormLoss_Struc(1, im) = trapz(x_PIDR, NormLoss_Struc_edpj); 
    end
   
    % =====================================================================
    % NON-STRUCTURAL DRIFT-SENSITIVE LOSS ASSESSMENT - HAZUS APPROACH
    % =====================================================================
    % HAZUS treats drift-sensitive non-structural components as a single
    % category or assembly
    
    % Select fragility parameters based on building height
    if Num_stories == 1
        medianEDP_NSDrift = mediandemandNonstructuralDrift(1, :) / 216; 
    elseif Num_stories == 5
        medianEDP_NSDrift = mediandemandNonstructuralDrift(2, :) / 540; 
    elseif Num_stories > 5
        medianEDP_NSDrift = mediandemandNonstructuralDrift(3, :) / 1123;
    end  
    
    dispEDP_NSDrift = disp_NonstructuralDrift;
    CostData_NSDrift = CostData(3, :);  % Cost ratios for drift-sensitive components

    % Loop through intensity measure levels
    for im = 1:length(SaCalc)

        % Use same demand model as structural (drift-based)
        MD_PIDR = PIDR_DemandPar(2) * log(SaCalc(im)) + PIDR_DemandPar(1);
        PDF_PIDR_IM = lognpdf(x_PIDR, MD_PIDR, SIGMA);

        % Initialize damage state calculation arrays
        Frag_DS_NSD = zeros(length(x_PIDR), length(medianEDP_NSDrift));
        Prob_ds_NSD = zeros(length(x_PIDR), length(medianEDP_NSDrift));
        Loss_EDP_NSD = zeros(length(x_PIDR), length(medianEDP_NSDrift));        

        % Calculate fragility functions for each damage state
        for ds = 1:length(medianEDP_NSDrift)
            Frag_DS_NSD(:, ds) = logncdf(x_PIDR, log(medianEDP_NSDrift(1, ds)), dispEDP_NSDrift);                
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
        
        % Sum losses across damage states and integrate
        NormLoss_NonStrucDrift_edpj = sum(Loss_EDP_NSD, 2);
        NormLoss_NonStrucDrift(1, im) = trapz(x_PIDR, NormLoss_NonStrucDrift_edpj); 
    end
       
    % =====================================================================
    % NON-STRUCTURAL ACCELERATION-SENSITIVE LOSS ASSESSMENT - HAZUS APPROACH
    % =====================================================================
    % HAZUS treats acceleration-sensitive components as a single category
    % or assembly
    
    % Select fragility parameters based on building height
    if Num_stories == 1
        medianEDP_NSAcc = mediandemandNonstructuralAcc(1, :); 
    elseif Num_stories == 5
        medianEDP_NSAcc = mediandemandNonstructuralAcc(2, :); 
    elseif Num_stories > 5
        medianEDP_NSAcc = mediandemandNonstructuralAcc(3, :);
    end  
    
    dispEDP_NSAcc = disp_NonstructuralAcc;
    CostData_NSAcc = CostData(2, :);    % Cost ratios for acceleration-sensitive components

    % Develop demand model: Sa vs Maximum Peak Floor Acceleration
    md_PFA = fitlm(log(IM_GM(bd, :)), log(Max_PFA(1, :)));
    PFA_DemandPar(1:2) = md_PFA.Coefficients{1:2, 1};      % [intercept, slope]
    PFA_DemandPar(3) = md_PFA.RMSE;                         % Residual standard error
    SIGMA_PFA = PFA_DemandPar(3);                           % Log-standard deviation

    % Loop through intensity measure levels
    for im = 1:length(SaCalc)

        % Calculate mean log demand for current intensity level
        MD_PFA = PFA_DemandPar(2) * log(SaCalc(im)) + PFA_DemandPar(1);
        
        % Calculate probability density function of acceleration given Sa
        PDF_PFA_IM = lognpdf(x_PFA, MD_PFA, SIGMA_PFA);

        % Initialize damage state calculation arrays
        Frag_DS_NSA = zeros(length(x_PFA), length(medianEDP_NSAcc));
        Prob_ds_NSA = zeros(length(x_PFA), length(medianEDP_NSAcc));
        Loss_EDP_NSA = zeros(length(x_PFA), length(medianEDP_NSAcc));            
            
        % Calculate fragility functions for each damage state
        for ds = 1:length(medianEDP_NSAcc)                
            Frag_DS_NSA(:, ds) = logncdf(x_PFA, log(medianEDP_NSAcc(1, ds)), dispEDP_NSAcc);                    
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
        
        % Sum losses across damage states and integrate
        NormLoss_NonStrucAcc_edpj = sum(Loss_EDP_NSA, 2); 
        NormLoss_NonStrucAcc(1, im) = trapz(x_PFA, NormLoss_NonStrucAcc_edpj); 
    end
    
    %% Normalize Loss Values
    % Convert from percentage to decimal form (HAZUS uses percentage values)
    NormLoss_Struc = NormLoss_Struc(1, :) ./ 100;
    NormLoss_NonStrucDrift = NormLoss_NonStrucDrift(1, :) ./ 100;
    NormLoss_NonStrucAcc = NormLoss_NonStrucAcc(1, :) ./ 100;
    
    % Calculate total normalized loss (sum of all component categories)
    NormTotLoss = NormLoss_Struc + NormLoss_NonStrucDrift + NormLoss_NonStrucAcc;

    % =====================================================================
    % DEMOLITION LOSS ASSESSMENT
    % =====================================================================
    % Same approach as FEMA P-58: based on residual drift thresholds
    
    % Define residual drift range for PDF calculation
    x_RSDR_pdf = [0.0002:0.0002:0.06]';
    
    % Demolition decision fragility curve from Ramirez and Miranda (2012)
    theta_hat_Demolish = 0.015;     % Median residual drift capacity (1.5%)
    beta_hat_Demolish = 0.3;        % Dispersion (log-standard deviation)
    
    % Calculate probability of demolition given residual drift
    PD_RSDR = normcdf((log(x_RSDR_pdf/theta_hat_Demolish))/beta_hat_Demolish);
    
    % Initialize probability of demolition given intensity measure
    PD_IM_NC = zeros(length(SaCalc), 1);  

    % Loop through each intensity measure level
    for im = 1:length(SaCalc)
    
        % Develop demand model for maximum residual drift
        md_RSDR = fitlm(log(IM_GM(bd, :)), log(Max_RSDR(1, :)));
        RSDR_DemandPar(1:2) = md_RSDR.Coefficients{1:2, 1};
        RSDR_DemandPar(3) = md_RSDR.RMSE;
        
        SIGMA_RSDR = RSDR_DemandPar(3);
        MD_RSDR = RSDR_DemandPar(2) * log(SaCalc(im)) + RSDR_DemandPar(1);    
    
        % Calculate PDF of residual drift given intensity measure
        dP_RIDR_IM = lognpdf(x_RSDR_pdf, MD_RSDR, SIGMA_RSDR);
    
        % Combine demolition probability with drift demand PDF
        PD_RSDR_x_dP_RSDR_IM = PD_RSDR(:, 1) .* dP_RIDR_IM(:, 1);    
    
        % Integrate to get probability of demolition given IM
        PD_IM_NC(im, 1) = trapz(x_RSDR_pdf, PD_RSDR_x_dP_RSDR_IM);
    end
    
    % Transpose for consistency
    PD_IM_NC = PD_IM_NC';
    
    % =====================================================================
    % REPLACEMENT COST CALCULATION
    % =====================================================================
    % Calculate building replacement cost
    
    Num_bays = 5;                       % Number of bays in each direction
    Floor_area = (Bay_width * Num_bays) * (Bay_width * Num_bays); % Floor area [ft²]
    Per_sqft_Cost = 250;                % Unit cost [$/ft²] from Hwang and Lignos (2017)
    Tot_Replac_Cost = Num_stories * Floor_area * Per_sqft_Cost;   % Total replacement cost [$]
    
    % Demolition cost (10% of replacement cost)
    DemoLoss = 0.1 * Tot_Replac_Cost;
    Tot_Loss_Demo = Tot_Replac_Cost + DemoLoss;
    
    % =====================================================================
    % EXPECTED LOSS CALCULATION WITH DEMOLITION CONSIDERATION
    % =====================================================================
    % Apply demolition probability to calculate final expected losses
    
    % Structural repair loss (conditional on no demolition)
    Norm_Loss_Struc(1, :) = NormLoss_Struc(1, :) .* (1. - PD_IM_NC(1, :));
    Loss_Struc(1, :) = Norm_Loss_Struc(1, :) .* Tot_Replac_Cost;
        
    % Non-structural drift-sensitive repair loss (conditional on no demolition)
    Norm_Loss_NonStruc_Drift(1, :) = NormLoss_NonStrucDrift(1, :) .* (1. - PD_IM_NC(1, :));
    Loss_NonStruc_Drift(1, :) = Norm_Loss_NonStruc_Drift(1, :) .* Tot_Replac_Cost;
    
    % Non-structural acceleration-sensitive repair loss (conditional on no demolition)
    Norm_Loss_NonStruc_Acc(1, :) = NormLoss_NonStrucAcc(1, :) .* (1. - PD_IM_NC(1, :));
    Loss_NonStruc_Acc(1, :) = Norm_Loss_NonStruc_Acc(1, :) .* Tot_Replac_Cost;
    
    % Demolition loss (conditional on demolition decision)
    Loss_Demo(1, :) = Tot_Loss_Demo * (PD_IM_NC(1, :)); 
    Norm_Loss_Demo(1, :) = Loss_Demo(1, :) ./ Tot_Replac_Cost;
    
    % Total expected loss (repair + demolition)
    Norm_Tot_Loss(1, :) = NormTotLoss(1, :) .* (1. - PD_IM_NC(1, :)) + Norm_Loss_Demo;  
    Loss_Tot(1, :) = Norm_Tot_Loss(1, :) .* Tot_Replac_Cost;

    % =====================================================================
    % SEISMIC HAZARD INPUT AND PROCESSING
    % =====================================================================
    % Load USGS seismic hazard data and interpolate for building-specific periods
    
    Hazard_file = 'USGSHazard_data.csv';
    HazardData = readmatrix(Hazard_file);
    AllLamda_IM = HazardData(:, 2:12);   % Exceedance rates for different periods
    IM_values = HazardData(:, 1);        % Spectral acceleration values [g]
    
    % Extract hazard curves for different structural periods
    a = AllLamda_IM(:, 2);   % T = 0.1s
    b = AllLamda_IM(:, 3);   % T = 0.2s
    c = AllLamda_IM(:, 4);   % T = 0.3s
    d = AllLamda_IM(:, 5);   % T = 0.5s
    e = AllLamda_IM(:, 6);   % T = 0.75s
    f = AllLamda_IM(:, 7);   % T = 1.0s
    g = AllLamda_IM(:, 8);   % T = 2.0s
    h = AllLamda_IM(:, 9);   % T = 3.0s
    i = AllLamda_IM(:, 10);  % T = 4.0s
    j = AllLamda_IM(:, 11);  % T = 5.0s
    
    % Initialize hazard curve for target period
    Lamda_IM_target = zeros(length(IM_values), length(Target_time_period));
    
    % Interpolate hazard curve for building's fundamental period
    % Linear interpolation in log space between adjacent period values
    if Target_time_period(bd, 1) >= 0.1 && Target_time_period(bd, 1) < 0.2
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([0.1, 0.2], [log(a(nn, 1)), log(b(nn, 1))], Target_time_period(bd, 1), 'linear');
        end

    elseif Target_time_period(bd, 1) >= 0.2 && Target_time_period(bd, 1) < 0.3
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([0.2, 0.3], [log(b(nn, 1)), log(c(nn, 1))], Target_time_period(bd, 1), 'linear');
        end

    elseif Target_time_period(bd, 1) >= 0.3 && Target_time_period(bd, 1) < 0.5
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([0.3, 0.5], [log(c(nn, 1)), log(d(nn, 1))], Target_time_period(bd, 1), 'linear');
        end

    elseif Target_time_period(bd, 1) >= 0.5 && Target_time_period(bd, 1) < 0.75
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([0.5, 0.75], [log(d(nn, 1)), log(e(nn, 1))], Target_time_period(bd, 1), 'linear');
        end

    elseif Target_time_period(bd, 1) >= 0.75 && Target_time_period(bd, 1) < 1
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([0.75, 1], [log(e(nn, 1)), log(f(nn, 1))], Target_time_period(bd, 1), 'linear');
        end

    elseif Target_time_period(bd, 1) >= 1 && Target_time_period(bd, 1) < 2
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([1, 2], [log(f(nn, 1)), log(g(nn, 1))], Target_time_period(bd, 1), 'linear');
        end

    elseif Target_time_period(bd, 1) >= 2 && Target_time_period(bd, 1) < 3
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([2, 3], [log(g(nn, 1)), log(h(nn, 1))], Target_time_period(bd, 1), 'linear');
        end

    elseif Target_time_period(bd, 1) >= 3 && Target_time_period(bd, 1) < 4
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([3, 4], [log(h(nn, 1)), log(i(nn, 1))], Target_time_period(bd, 1), 'linear');
        end

    elseif Target_time_period(bd, 1) >= 4 && Target_time_period(bd, 1) < 5
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn, bd) = interp1([4, 5], [log(i(nn, 1)), log(j(nn, 1))], Target_time_period(bd, 1), 'linear');
        end
    end
    
    % Format hazard curve as [Sa, λ(Sa)] pairs
    hazardCurveSa = [IM_values, exp(Lamda_IM_target(:, bd))];
    
    % Calculate derivative of hazard curve (dλ/dSa) using numerical differentiation
    PDF_temp = zeros(length(SaCalc), 1);
    hh = 0.0000001;  % Small increment for numerical differentiation
    
    for kk = 1:length(SaCalc)
        % Calculate derivative using two-point formula
        % f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
        
        % Interpolate λ(Sa + hh)
        fxp = interp1(hazardCurveSa(:, 1), hazardCurveSa(:, 2), SaCalc(kk) + hh, 'spline'); 
        % Interpolate λ(Sa - hh)
        fxn = interp1(hazardCurveSa(:, 1), hazardCurveSa(:, 2), SaCalc(kk) - hh, 'spline');  
        % Calculate derivative
        PDF_temp(kk) = abs((fxp - fxn)) / (2 * hh);             
    end                                      
    
    % Smooth the derivative to reduce numerical noise
    dlambdaSa(:, bd) = smooth(PDF_temp);

    % =====================================================================
    % EXPECTED ANNUAL LOSS (EAL) COMPUTATION
    % =====================================================================
    % Integrate loss functions with hazard curve derivatives to calculate EAL
    % EAL = ∫ Loss(Sa) × |dλ/dSa| dSa
    
    % Calculate convolution of losses with hazard curve derivatives
    Conv_loss_Struc = Loss_Struc(1, :)' .* dlambdaSa(:, bd);
    Conv_loss_NonStrucDrift = Loss_NonStruc_Drift(1, :)' .* dlambdaSa(:, bd);
    Conv_loss_NonStrucAcc = Loss_NonStruc_Acc(1, :)' .* dlambdaSa(:, bd);
    Conv_loss_Tot = Loss_Tot(1, :)' .* dlambdaSa(:, bd);
    Conv_loss_Demo = Loss_Demo(1, :)' .* dlambdaSa(:, bd);

    % Calculate convolution of normalized losses
    Conv_Norm_loss_Struc = Norm_Loss_Struc(1, :)' .* dlambdaSa(:, bd);
    Conv_Norm_loss_NonStrucDrift = Norm_Loss_NonStruc_Drift(1, :)' .* dlambdaSa(:, bd);
    Conv_Norm_loss_NonStrucAcc = Norm_Loss_NonStruc_Acc(1, :)' .* dlambdaSa(:, bd);
    Conv_Norm_loss_Tot = Norm_Tot_Loss(1, :)' .* dlambdaSa(:, bd);
    Conv_Norm_loss_Demo = Norm_Loss_Demo(1, :)' .* dlambdaSa(:, bd);

    % Integrate to calculate Expected Annual Loss (EAL) values [$]
    EAL_Tot(bd, 1) = trapz(SaCalc, Conv_loss_Tot);
    EAL_Demo(bd, 1) = trapz(SaCalc, Conv_loss_Demo);
    EAL_Struc(bd, 1) = trapz(SaCalc, Conv_loss_Struc);
    EAL_NonStrucDrift(bd, 1) = trapz(SaCalc, Conv_loss_NonStrucDrift);
    EAL_NonStrucAcc(bd, 1) = trapz(SaCalc, Conv_loss_NonStrucAcc);

    % Calculate Normalized EAL (fraction of replacement cost) [-]
    Norm_EAL_Tot(bd, 1) = trapz(SaCalc, Conv_Norm_loss_Tot);
    Norm_EAL_Demo(bd, 1) = trapz(SaCalc, Conv_Norm_loss_Demo);
    Norm_EAL_Struc(bd, 1) = trapz(SaCalc, Conv_Norm_loss_Struc);
    Norm_EAL_NonStrucDrift(bd, 1) = trapz(SaCalc, Conv_Norm_loss_NonStrucDrift);
    Norm_EAL_NonStrucAcc(bd, 1) = trapz(SaCalc, Conv_Norm_loss_NonStrucAcc);
end

% =========================================================================
% OUTPUT FORMATTING AND EXPORT
% =========================================================================
% Compile all results into a comprehensive output matrix

% Create output matrix with all calculated values
All_Output(:, 1) = ID(:, 1);                           % Building ID
All_Output(:, 2) = EAL_Struc(:, 1);                    % Structural EAL [$]
All_Output(:, 3) = EAL_NonStrucDrift(:, 1);            % Non-structural drift EAL [$]
All_Output(:, 4) = EAL_NonStrucAcc(:, 1);              % Non-structural acceleration EAL [$]
All_Output(:, 5) = EAL_Demo(:, 1);                     % Demolition EAL [$]
All_Output(:, 6) = EAL_Tot(:, 1);                      % Total EAL [$]
All_Output(:, 7) = Norm_EAL_Struc(:, 1);               % Normalized structural EAL [-]
All_Output(:, 8) = Norm_EAL_NonStrucDrift(:, 1);       % Normalized non-structural drift EAL [-]
All_Output(:, 9) = Norm_EAL_NonStrucAcc(:, 1);         % Normalized non-structural acceleration EAL [-]
All_Output(:, 10) = Norm_EAL_Demo(:, 1);               % Normalized demolition EAL [-]
All_Output(:, 11) = Norm_EAL_Tot(:, 1);                % Normalized total EAL [-]

% =========================================================================
% OPTIONAL: EXPORT TO EXCEL (COMMENTED OUT)
% =========================================================================
% Uncomment the following lines to save results to Excel file

% % Specify the Excel file name
% excelFileName = 'HAZUS_NormEAL_MAX_07-09-24.xlsx';
% 
% % Specify the sheet name
% sheetName = 'Sheet1';
% 
% % Write the data to Excel
% xlswrite(excelFileName, All_Output, sheetName);
% disp('Data written to Excel successfully.');

% Display comprehensive analysis summary
disp('=================================================================');
disp('                   HAZUS ANALYSIS COMPLETE                      ');
disp('=================================================================');

disp('=================================================================');
disp('Results stored in variable: All_Output');
disp('Matrix columns:');
disp('  1: Building ID');
disp('  2: Structural EAL [$]');
disp('  3: Non-structural drift EAL [$]');
disp('  4: Non-structural acceleration EAL [$]');
disp('  5: Demolition EAL [$]');
disp('  6: Total EAL [$]');
disp('  7: Normalized structural EAL [-]');
disp('  8: Normalized non-structural drift EAL [-]');
disp('  9: Normalized non-structural acceleration EAL [-]');
disp(' 10: Normalized demolition EAL [-]');
disp(' 11: Normalized total EAL [-]');
disp('=================================================================');

% =========================================================================
% END OF HAZUS SEISMIC LOSS ASSESSMENT
% =========================================================================
