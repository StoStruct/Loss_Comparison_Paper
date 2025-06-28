%% FEMA P-58 Loss Calculation Main Code
% =========================================================================
% PROBABILISTIC SEISMIC LOSS ASSESSMENT FOR STEEL MOMENT-RESISTING FRAMES
% =========================================================================
%
% Description:
%   This script implements the FEMA P-58 methodology for probabilistic
%   seismic loss assessment of steel moment-resisting frame buildings.
%   The assessment includes:
%   - Structural components (columns, beams, connections)
%   - Non-structural drift-sensitive components (partitions, facades)
%   - Non-structural acceleration-sensitive components (equipment, contents)
%   - Demolition losses based on residual drift thresholds
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
%   FEMA P-58-1 (2018). Seismic Performance Assessment of Buildings,
%   Volume 1 - Methodology. Federal Emergency Management Agency.
%
% Input Files Required:
%   - GuanDataBase-IMs.csv: Intensity measures (Sa) for 240 ground motions
%   - Building geometry files in Guan_database/BuildingDesigns/
%   - Structural response files (PSDR, PFA, RSDR) in Guan_database/
%   - Fragility and cost function files in FragilityCostFunctions/
%   - USGSHazard_data.csv: Seismic hazard curves for different periods
%
% Output:
%   - Expected Annual Loss (EAL) values for each building [$]
%   - Normalized EAL values (as fraction of replacement cost) [-]
%   - Component-level breakdown of losses
%
% Usage:
%   Run this script directly. Ensure all required data files are in the
%   correct directory structure as specified in the file paths below.

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
SaCalc = [0.001, 0.01, 0.05, 0.1:0.01:2]'; % Sa values from 0.001g to 2.0g

%% Initialize Output Arrays
% Pre-allocate arrays for computational efficiency
dlambdaSa = zeros(length(SaCalc), length(ID));     % Hazard curve derivatives

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

    % Load response data
    PSDR = readmatrix(filePSDR);         % Peak story drift ratios (dimensionless)
    PFA = readmatrix(filePFA);           % Peak floor accelerations (in/s²)

    % Convert PFA from in/s² to g (gravity units)
    PFA = PFA * 0.0254 * 0.10197;       % Conversion factor: in/s² to g

    RSDR = readmatrix(fileRSDR);         % Residual story drift ratios (dimensionless)
    PSDR_trans = PSDR';                  % Transpose for processing
    Geometry = readmatrix(fileGeometry); % Building geometry parameters

    % =====================================================================
    % STRUCTURAL COMPONENT LOSS ASSESSMENT
    % =====================================================================
    % Load fragility and cost functions for structural components
    fileFragCostComp_Structure = fullfile('FragilityCostFunctions', ...
        'StructuralComponents', 'FragilityCostFunction_StructuralComponents.csv');
    FragCost_data_Structure = readmatrix(fileFragCostComp_Structure);

    % Extract building configuration parameters
    Num_LFRS_Bays = Geometry(1, 2);     % Number of lateral force resisting system bays
    Num_stories = Geometry(1, 1);       % Number of stories
    Bay_width = Geometry(1, 6);         % Bay width [ft]
    Num_bays = 5;                       % Total number of bays (5x5 grid)
    Floor_area = (Bay_width * Num_bays)^2; % Floor area [ft²]

    % Calculate component quantities based on LFRS configuration
    % Different configurations have different numbers of moment connections
    switch Num_LFRS_Bays
        case 1  % Single bay LFRS
            SMRF_Cols = 8;               % Special moment-resisting frame columns
            RBS_1sided = 8;              % Reduced beam section connections (1-sided)
            RBS_2sided = 0;              % Reduced beam section connections (2-sided)
        case 3  % Three bay LFRS
            SMRF_Cols = 16;
            RBS_1sided = 8;
            RBS_2sided = 8;
        case 5  % Five bay LFRS
            SMRF_Cols = 24;
            RBS_1sided = 0;
            RBS_2sided = 20;
    end

    % Calculate total structural components
    TotalColumns = 36;                   % Total columns in 6x6 grid
    TotalBeams = 25;                     % Total beams in 5x5 bays
    SMRF_Beams = Num_LFRS_Bays * 4;     % Moment frame beams
    Num_GravityBeams = TotalBeams - SMRF_Beams; % Gravity-only beams
    Num_ShearTabs = Num_GravityBeams * 2; % Shear tab connections

    % Component configuration matrix: [quantity, floor_pattern]
    % floor_pattern: 1=first floor only, 2=every 2nd floor, 3=all floors
    component_config = [
        SMRF_Cols,      1;  % Component 1: SMRF Columns (For base plate connections)
        TotalColumns,   2;  % Component 2: All Columns (For splices)
        Num_ShearTabs,  3;  % Component 3: Shear Tab connections
        RBS_1sided,     3;  % Component 4: RBS 1-sided connections
        RBS_2sided,     3;  % Component 5: RBS 2-sided connections
        Floor_area,     3   % Component 6: Floor System (per unit area)
        ];

    % Initialize structural loss arrays
    AllCompLoss_Struc = zeros(length(FragCost_data_Structure), 1);
    TotCompLoss_Struc = zeros(1, length(SaCalc));

    % Define drift range for integration
    x_PIDR = (0.0002:0.0002:3)';        % Peak inter-story drift range

    % Loop through each intensity measure level
    for im = 1:length(SaCalc)

        % Loop through each structural component type
        for comp = 1:6  % 6 structural component types

            % Extract fragility data from the structured matrix
            % Data is organized in sets of 4 columns per component
            data_col = (comp-1) * 4 + 1;
            medianEDP = FragCost_data_Structure(:, data_col);     % Median drift capacity
            dispEDP = FragCost_data_Structure(:, data_col + 1);   % Dispersion (log-std)
            medianCOST = FragCost_data_Structure(:, data_col + 2); % Median repair cost

            % Filter out invalid damage states (NaN or zero values)
            valid_idx = ~isnan(medianEDP) & medianEDP > 0;
            medianEDP = medianEDP(valid_idx);
            dispEDP = dispEDP(valid_idx);
            medianCOST = medianCOST(valid_idx);

            % Skip component if no valid damage states
            if isempty(medianEDP)
                continue;
            end

            % Get component quantity and floor distribution pattern
            quantity = component_config(comp, 1);
            floor_pattern = component_config(comp, 2);

            % Special case: skip all columns for single-story buildings
            if comp == 2 && Num_stories == 1
                AllCompLoss_Struc(data_col, 1) = 0;
                continue;
            end

            % Determine which floors to analyze based on component type
            switch floor_pattern
                case 1  % First floor only (base connections)
                    floors = 1;
                case 2  % Every 2nd floor (intermediate)
                    floors = 1:2:height(PSDR);
                case 3  % All floors (distributed components)
                    floors = 1:height(PSDR);
            end

            % Initialize floor-level loss array
            SumLoss_EDP_FLR_struct = zeros(height(PSDR), 1);

            % Calculate losses for each relevant floor
            for flr = floors
                % Develop demand model: log(PSDR) vs log(Sa) regression
                md_PIDR = fitlm(log(IM_GM(bd, :)), log(PSDR(flr, :)));
                PIDR_DemandPar(1:2) = md_PIDR.Coefficients{1:2, 1}; % [intercept, slope]
                PIDR_DemandPar(3) = md_PIDR.RMSE;                   % Residual standard error

                % Extract demand model parameters
                SIGMA = PIDR_DemandPar(3);                          % Log-standard deviation
                MD_PIDR = PIDR_DemandPar(2) * log(SaCalc(im)) + PIDR_DemandPar(1); % Mean log(PSDR)

                % Calculate probability density function of PSDR given Sa
                PDF_PIDR_IM = lognpdf(x_PIDR, MD_PIDR, SIGMA);

                % Initialize arrays for damage state calculations
                Loss_edpj_struct = zeros(length(x_PIDR), 1);
                Frag_DS = zeros(length(x_PIDR), length(medianEDP));
                Prob_ds = zeros(length(x_PIDR), length(medianEDP));
                Loss_EDP = zeros(length(x_PIDR), length(medianEDP));

                % Calculate fragility functions for each damage state
                for ds = 1:length(medianEDP)
                    % Lognormal CDF for exceedance probability
                    Frag_DS(:, ds) = logncdf(x_PIDR, log(medianEDP(ds)), dispEDP(ds));
                end

                % Calculate probability of being in each damage state
                % P(DS=i) = F(DS≥i) - F(DS≥i+1)
                if length(medianEDP) == 1
                    Prob_ds(:, 1) = Frag_DS(:, 1);
                else
                    for p = 1:length(medianEDP)-1
                        Prob_ds(:, p) = Frag_DS(:, p) - Frag_DS(:, p+1);
                    end
                    Prob_ds(:, end) = Frag_DS(:, end);
                end

                % Calculate expected loss for each damage state
                for p = 1:length(medianEDP)
                    Loss_EDP(:, p) = medianCOST(p) * quantity * (Prob_ds(:, p) .* PDF_PIDR_IM);
                end

                % Sum losses across all damage states
                Loss_edpj_struct(:, 1) = sum(Loss_EDP, 2);

                % Integrate over the demand range to get expected loss
                SumLoss_EDP_FLR_struct(flr, 1) = trapz(x_PIDR, Loss_edpj_struct);
            end

            % Store total component loss (sum across all floors)
            AllCompLoss_Struc(data_col, 1) = sum(SumLoss_EDP_FLR_struct(:, 1));
        end

        % Calculate total structural loss for this intensity level
        TotCompLoss_Struc(1, im) = sum(AllCompLoss_Struc(:, 1));
    end

    % =====================================================================
    % NON-STRUCTURAL DRIFT-SENSITIVE COMPONENT LOSS ASSESSMENT
    % =====================================================================
    % Load fragility and cost data for drift-sensitive non-structural components
    fileMedianDemand_NonStructureDrift = strcat(['FragilityCostFunctions', ...
        '\NonStructuralDriftSensitiveComponents'], '\', ...
        'MedianDemand_NonStructuralDriftSensitive.csv');

    fileDispersionDemand_NonStructuralDrift = strcat(['FragilityCostFunctions', ...
        '\NonStructuralDriftSensitiveComponents'], '\', ...
        'DispersionDemand_NonStructuralDriftSensitive.csv');

    fileMedianCost_NonStructureDrift = strcat(['FragilityCostFunctions', ...
        '\NonStructuralDriftSensitiveComponents'], '\', ...
        'MedianCost_NonStructuralDriftSensitive.csv');

    fileQuantities_NonStructureDrift = strcat(['FragilityCostFunctions', ...
        '\NonStructuralDriftSensitiveComponents'], '\', ...
        'Quantity_NonStructuralDriftSensitive - Copy.csv');

    % Load non-structural drift-sensitive component data
    mediandemand = readmatrix(fileMedianDemand_NonStructureDrift);
    dispdemand = readmatrix(fileDispersionDemand_NonStructuralDrift);
    mediancost = readmatrix(fileMedianCost_NonStructureDrift);
    quantities = readmatrix(fileQuantities_NonStructureDrift);

    % Component types: partitions, facades, stairs, elevator systems etc
    Num_NonStructDriftComp = 4;
    AllCompLoss_NonStrucDrift = zeros(Num_NonStructDriftComp, 1);
    TotCompLoss_NonStrucDrift = zeros(1, length(SaCalc));

    % Define drift range for integration (same as structural)
    x_PIDR = (0.0002:0.0002:3)';

    % Loop through intensity measure levels
    for im = 1:length(SaCalc)

        % Loop through each non-structural drift-sensitive component
        for comp = 1:Num_NonStructDriftComp

            % Define number of damage states for each component type
            % Different components have different vulnerability characteristics
            if comp == 1
                num_ds = 2;     % Component 1: 2 damage states
            elseif comp == 2
                num_ds = 3;     % Component 2: 3 damage states
            elseif comp == 3
                num_ds = 1;     % Component 3: 1 damage state
            elseif comp == 4
                num_ds = 3;     % Component 4: 3 damage states
            end

            % Extract fragility parameters for current component
            medianEDP = (mediandemand(1:num_ds, comp))';   % Median drift capacities
            dispEDP = (dispdemand(1:num_ds, comp))';       % Dispersions

            % Assign component quantities and costs based on bay width
            % Buildings with different bay widths have different quantities
            Bay_width = Geometry(1, 6);

            if Bay_width == 20
                CompQuantity = quantities(1, comp);
                medianCOST = (mediancost(1:num_ds, comp))';
            elseif Bay_width == 30
                CompQuantity = quantities(2, comp);
                medianCOST = (mediancost(1:num_ds, comp+Num_NonStructDriftComp))';
            elseif Bay_width == 40
                CompQuantity = quantities(3, comp);
                medianCOST = (mediancost(1:num_ds, comp+Num_NonStructDriftComp*2))';
            end

            % Initialize floor-level loss array
            SumLoss_EDP_FLR_NSD = zeros(height(PSDR), 1);

            % Calculate loss for each floor
            for flr = 1:height(PSDR)

                % Develop demand model for current floor
                md_PIDR = fitlm(log(IM_GM(bd, :)), log(PSDR(flr, :)));
                PIDR_DemandPar(1:2) = md_PIDR.Coefficients{1:2, 1};
                PIDR_DemandPar(3) = md_PIDR.RMSE;

                SIGMA = PIDR_DemandPar(3);
                MD_PIDR = PIDR_DemandPar(2) * log(SaCalc(im)) + PIDR_DemandPar(1);

                % Calculate PDF of drift given intensity measure
                PDF_PIDR_IM_NSD = lognpdf(x_PIDR, MD_PIDR, SIGMA);

                % Initialize damage state calculation arrays
                Loss_edpj_NSD = zeros(length(x_PIDR), 1);
                Frag_DS_NSD = zeros(length(x_PIDR), num_ds);
                Prob_ds_NSD = zeros(length(x_PIDR), num_ds);
                Loss_EDP_NSD = zeros(length(x_PIDR), num_ds);

                % Calculate fragility functions for each damage state
                for ds = 1:num_ds
                    Frag_DS_NSD(:, ds) = logncdf(x_PIDR, log(medianEDP(1, ds)), dispEDP(1, ds));
                end

                % Calculate probability of being in each damage state
                if num_ds == 1
                    Prob_ds_NSD(:, 1) = Frag_DS_NSD(:, 1);
                else
                    for p = 1:num_ds-1
                        Prob_ds_NSD(:, p) = Frag_DS_NSD(:, p) - Frag_DS_NSD(:, p+1);
                        Prob_ds_NSD(:, p+1) = Frag_DS_NSD(:, p+1);
                    end
                end

                % Calculate expected loss for each damage state
                for p = 1:num_ds
                    Loss_EDP_NSD(:, p) = medianCOST(1, p) * CompQuantity * (Prob_ds_NSD(:, p) .* PDF_PIDR_IM_NSD);
                end

                % Sum losses across damage states
                Loss_edpj_NSD(:, 1) = sum(Loss_EDP_NSD, 2);

                % Integrate to get floor-level expected loss
                SumLoss_EDP_FLR_NSD(flr, 1) = trapz(x_PIDR, Loss_edpj_NSD);
            end

            % Sum losses across all floors for this component
            AllCompLoss_NonStrucDrift(comp, 1) = sum(SumLoss_EDP_FLR_NSD(:, 1));
        end

        % Calculate total non-structural drift-sensitive loss
        TotCompLoss_NonStrucDrift(1, im) = sum(AllCompLoss_NonStrucDrift(:, 1));
    end

    % =====================================================================
    % NON-STRUCTURAL ACCELERATION-SENSITIVE COMPONENT LOSS ASSESSMENT
    % =====================================================================
    % These components include equipment, contents, and building systems
    % that are sensitive to floor accelerations 

    % ----------------------------- FLOOR-LEVEL COMPONENTS ---------------
    % Load fragility and cost data for acceleration-sensitive components
    fileMedianDemand_NonStructureAcc = strcat(['FragilityCostFunctions', ...
        '\NonStructuralAccSensitiveComponents'], '\', ...
        'MedianDemand_NonStructuralAccSensitive.csv');

    fileDispersionDemand_NonStructuralAcc = strcat(['FragilityCostFunctions', ...
        '\NonStructuralAccSensitiveComponents'], '\', ...
        'DispersionDemand_NonStructuralAccSensitive.csv');

    fileMedianCost_NonStructureAcc = strcat(['FragilityCostFunctions', ...
        '\NonStructuralAccSensitiveComponents'], '\', ...
        'MedianCost_NonStructuralAccSensitive.csv');

    fileQuantities_NonStructureAcc = strcat(['FragilityCostFunctions', ...
        '\NonStructuralAccSensitiveComponents'], '\', ...
        'Quantity_NonStructuralAccSensitive - Copy.csv');

    % Load acceleration-sensitive component data
    mediandemand = readmatrix(fileMedianDemand_NonStructureAcc, 'EmptyValue', 0);
    dispdemand = readmatrix(fileDispersionDemand_NonStructuralAcc, 'EmptyValue', 0);
    mediancost = readmatrix(fileMedianCost_NonStructureAcc, 'EmptyValue', 0);
    quantities = readmatrix(fileQuantities_NonStructureAcc);

    % 12 different types of acceleration-sensitive components
    Num_NonStructAccComp = 12;
    AllCompLoss_NonStrucAcc_FlrLvl = zeros(Num_NonStructAccComp, 1);
    TotCompLoss_NonStrucAcc_FlrLvl = zeros(1, length(SaCalc));

    % Define acceleration range for integration [g]
    x_PFA = (0.02:0.02:10)';

    % Loop through intensity measure levels
    for im = 1:length(SaCalc)

        % Loop through each acceleration-sensitive component type
        for comp = 1:Num_NonStructAccComp

            % Define number of damage states for each component type
            % Different equipment types have different vulnerability patterns
            if comp == 1
                num_ds = 2;
            elseif comp == 2
                num_ds = 1;
            elseif comp == 3
                num_ds = 3;
            elseif comp == 4
                num_ds = 2;
            elseif comp == 5
                num_ds = 2;
            elseif comp == 6
                num_ds = 2;
            elseif comp == 7
                num_ds = 1;
            elseif comp == 8
                num_ds = 1;
            elseif comp == 9
                num_ds = 1;
            elseif comp == 10
                num_ds = 2;
            elseif comp == 11
                num_ds = 2;
            elseif comp == 12
                num_ds = 1;
            end

            % Extract fragility parameters
            medianEDP = mediandemand(1:num_ds, comp);   % Median acceleration capacities [g]
            dispEDP = dispdemand(1:num_ds, comp);       % Dispersions

            % Assign quantities and costs based on bay width
            Bay_width = Geometry(1, 6);

            if Bay_width == 20
                CompQuantity = quantities(1, comp);
                medianCOST = mediancost(1:num_ds, comp);
            elseif Bay_width == 30
                CompQuantity = quantities(2, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccComp);
            elseif Bay_width == 40
                CompQuantity = quantities(3, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccComp*2);
            end

            % Initialize floor-level loss array
            SumLoss_EDP_FLR_NSA = zeros(height(PSDR), 1);

            % Calculate loss for each floor (excluding roof)
            for flr = 1:height(PFA)-1

                % Develop demand model: log(PFA) vs log(Sa) regression
                md_PFA = fitlm(log(IM_GM(bd, :)), log(PFA(flr, :)));
                PFA_DemandPar(1:2) = md_PFA.Coefficients{1:2, 1};
                PFA_DemandPar(3) = md_PFA.RMSE;

                SIGMA_PFA = PFA_DemandPar(3);
                MD_PFA = PFA_DemandPar(2) * log(SaCalc(im)) + PFA_DemandPar(1);

                % Calculate PDF of acceleration given intensity measure
                PDF_PFA_IM = lognpdf(x_PFA, MD_PFA, SIGMA_PFA);

                % Initialize damage state calculation arrays
                Loss_edpj_NSA = zeros(length(x_PFA), 1);
                Frag_DS_NSA = zeros(length(x_PFA), num_ds);
                Prob_ds_NSA = zeros(length(x_PFA), num_ds);
                Loss_EDP_NSA = zeros(length(x_PFA), num_ds);

                % Calculate fragility functions
                if num_ds == 1
                    Frag_DS_NSA(:, num_ds) = logncdf(x_PFA, log(medianEDP(num_ds, 1)), dispEDP(num_ds, 1));
                else
                    for ds = 1:num_ds
                        Frag_DS_NSA(:, ds) = logncdf(x_PFA, log(medianEDP(ds, 1)), dispEDP(ds, 1));
                    end
                end

                % Calculate probability of being in each damage state
                if num_ds == 1
                    Prob_ds_NSA(:, 1) = Frag_DS_NSA(:, 1);
                else
                    for p = 1:length(medianEDP)-1
                        Prob_ds_NSA(:, p) = Frag_DS_NSA(:, p) - Frag_DS_NSA(:, p+1);
                        Prob_ds_NSA(:, p+1) = Frag_DS_NSA(:, p+1);
                    end
                end

                % Calculate expected loss for each damage state
                for p = 1:length(medianEDP)
                    Loss_EDP_NSA(:, p) = medianCOST(p, 1) * CompQuantity * (Prob_ds_NSA(:, p) .* PDF_PFA_IM);
                end

                % Sum losses across damage states
                Loss_edpj_NSA(:, 1) = sum(Loss_EDP_NSA, 2);

                % Integrate to get floor-level expected loss
                SumLoss_EDP_FLR_NSA(flr, 1) = trapz(x_PFA, Loss_edpj_NSA);
            end

            % Sum losses across all floors for this component
            AllCompLoss_NonStrucAcc_FlrLvl(comp, 1) = sum(SumLoss_EDP_FLR_NSA(:, 1));
        end

        % Calculate total floor-level acceleration-sensitive loss
        TotCompLoss_NonStrucAcc_FlrLvl(1, im) = sum(AllCompLoss_NonStrucAcc_FlrLvl(:, 1));
    end

    % ----------------------------- BUILDING-LEVEL COMPONENTS -------------
    % These are components that respond to building-level accelerations
    % rather than individual floor accelerations

    % Load building-level acceleration-sensitive component data
    fileMedianDemand_NonStructureAcc_BuildingLevel = strcat('FragilityCostFunctions\NonStructuralAccSensitiveComponents\BuildingLevel', '\', 'MedianDemand_NonStructuralAccSensitive_Buildinglevel.csv');
    fileDispersionDemand_NonStructuralAcc_BuildingLevel = strcat('FragilityCostFunctions\NonStructuralAccSensitiveComponents\BuildingLevel', '\', 'DispersionDemand_NonStructuralAccSensitive_Buildinglevel.csv');
    fileMedianCost_NonStructureAcc_BuildingLevel = strcat('FragilityCostFunctions\NonStructuralAccSensitiveComponents\BuildingLevel', '\', 'MedianCost_NonStructuralAccSensitive_Buildinglevel.csv');
    fileQuantities_NonStructureAcc_BuildingLevel = strcat('FragilityCostFunctions\NonStructuralAccSensitiveComponents\BuildingLevel', '\', 'Quantity_NonStructuralAccSensitive_Buildinglevel - Copy.csv');

    % Load building-level component data
    mediandemand = readmatrix(fileMedianDemand_NonStructureAcc_BuildingLevel);
    dispdemand = readmatrix(fileDispersionDemand_NonStructuralAcc_BuildingLevel);
    mediancost = readmatrix(fileMedianCost_NonStructureAcc_BuildingLevel);
    quantities = readmatrix(fileQuantities_NonStructureAcc_BuildingLevel);

    % 5 types of building-level acceleration-sensitive components
    Num_NonStructAccCompBldg = 5;
    AllCompLoss_NonStrucAcc_BldgLvl = zeros(Num_NonStructAccCompBldg, 1);
    TotCompLoss_NonStrucAcc_BldgLvl = zeros(1, length(SaCalc));

    % Define acceleration range for building-level components [g]
    x_PFA = (0.02:0.02:5)';
    AllCompLoss_NonStrucAcc_BldgLvl_edpj = zeros(length(x_PFA), 1);

    % Loop through intensity measure levels
    for im = 1:length(SaCalc)

        % Loop through each building-level acceleration-sensitive component
        for comp = 1:Num_NonStructAccCompBldg

            % Define number of damage states for each building-level component
            if comp == 1
                num_ds = 4;     % Component 1: 4 damage states
            elseif comp == 2
                num_ds = 1;     % Component 2: 1 damage state
            elseif comp == 3
                num_ds = 1;     % Component 3: 1 damage state
            elseif comp == 4
                num_ds = 2;     % Component 4: 2 damage states
            elseif comp == 5
                num_ds = 1;     % Component 5: 1 damage state
            end

            % Extract fragility parameters
            medianEDP = mediandemand(1:num_ds, comp);   % Median acceleration capacities [g]
            dispEDP = dispdemand(1:num_ds, comp);       % Dispersions

            % Assign component quantities and costs based on building characteristics
            % Different combinations of stories and bay width require different quantities
            Bay_width = Geometry(1, 6);
            Num_Storey = Geometry(1, 1);

            % 15 different building configurations (5 story heights × 3 bay widths)
            if (Num_Storey == 1) && (Bay_width == 20)
                CompQuantity = quantities(1, comp);
                medianCOST = mediancost(1:num_ds, comp);
            elseif (Num_Storey == 1) && (Bay_width == 30)
                CompQuantity = quantities(2, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg);
            elseif (Num_Storey == 1) && (Bay_width == 40)
                CompQuantity = quantities(3, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*2);

            elseif (Num_Storey == 5) && (Bay_width == 20)
                CompQuantity = quantities(4, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*3);
            elseif (Num_Storey == 5) && (Bay_width == 30)
                CompQuantity = quantities(5, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*4);
            elseif (Num_Storey == 5) && (Bay_width == 40)
                CompQuantity = quantities(6, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*5);

            elseif (Num_Storey == 9) && (Bay_width == 20)
                CompQuantity = quantities(7, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*6);
            elseif (Num_Storey == 9) && (Bay_width == 30)
                CompQuantity = quantities(8, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*7);
            elseif (Num_Storey == 9) && (Bay_width == 40)
                CompQuantity = quantities(9, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*8);

            elseif (Num_Storey == 14) && (Bay_width == 20)
                CompQuantity = quantities(10, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*9);
            elseif (Num_Storey == 14) && (Bay_width == 30)
                CompQuantity = quantities(11, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*10);
            elseif (Num_Storey == 14) && (Bay_width == 40)
                CompQuantity = quantities(12, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*11);

            elseif (Num_Storey == 19) && (Bay_width == 20)
                CompQuantity = quantities(13, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*12);
            elseif (Num_Storey == 19) && (Bay_width == 30)
                CompQuantity = quantities(14, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*13);
            elseif (Num_Storey == 19) && (Bay_width == 40)
                CompQuantity = quantities(15, comp);
                medianCOST = mediancost(1:num_ds, comp+Num_NonStructAccCompBldg*14);
            end

            % Use first floor acceleration as representative building acceleration
            flr = 1;

            % Develop demand model for building-level acceleration
            md_PFA = fitlm(log(IM_GM(bd, :)), log(PFA(flr, :)));
            PFA_DemandPar(1:2) = md_PFA.Coefficients{1:2, 1};
            PFA_DemandPar(3) = md_PFA.RMSE;

            SIGMA_PFA = PFA_DemandPar(3);
            MD_PFA = PFA_DemandPar(2) * log(SaCalc(im)) + PFA_DemandPar(1);

            % Calculate PDF of acceleration given intensity measure
            PDF_PFA_IM_BldgLvl = lognpdf(x_PFA, MD_PFA, SIGMA_PFA);

            % Initialize damage state calculation arrays
            Loss_edpj_NSA_BldgLvl = zeros(length(x_PFA), 1);
            Frag_DS_NSA_BldgLvl = zeros(length(x_PFA), num_ds);
            Prob_ds_NSA_BldgLvl = zeros(length(x_PFA), num_ds);
            Loss_EDP_NSA_BldgLvl = zeros(length(x_PFA), num_ds);

            % Calculate fragility functions for each damage state
            for ds = 1:num_ds
                Frag_DS_NSA_BldgLvl(:, ds) = logncdf(x_PFA, log(medianEDP(ds, 1)), dispEDP(ds, 1));
            end

            % Calculate probability of being in each damage state
            if num_ds == 1
                Prob_ds_NSA_BldgLvl(:, 1) = Frag_DS_NSA_BldgLvl(:, 1);
            else
                for p = 1:length(medianEDP)-1
                    Prob_ds_NSA_BldgLvl(:, p) = Frag_DS_NSA_BldgLvl(:, p) - Frag_DS_NSA_BldgLvl(:, p+1);
                    Prob_ds_NSA_BldgLvl(:, p+1) = Frag_DS_NSA_BldgLvl(:, p+1);
                end
            end

            % Calculate expected loss for each damage state
            for p = 1:length(medianEDP)
                Loss_EDP_NSA_BldgLvl(:, p) = medianCOST(p, 1) * CompQuantity * (Prob_ds_NSA_BldgLvl(:, p) .* PDF_PFA_IM_BldgLvl);
            end

            % Sum losses across damage states
            Loss_edpj_NSA_BldgLvl(:, 1) = sum(Loss_EDP_NSA_BldgLvl, 2);

            % Integrate to get component-level expected loss
            AllCompLoss_NonStrucAcc_BldgLvl(comp, 1) = trapz(x_PFA, Loss_edpj_NSA_BldgLvl);
        end

        % Calculate total building-level acceleration-sensitive loss
        TotCompLoss_NonStrucAcc_BldgLvl(1, im) = sum(AllCompLoss_NonStrucAcc_BldgLvl(:, 1));
    end

    % Combine floor-level and building-level acceleration-sensitive losses
    TotCompLoss_NonStrucAcc = TotCompLoss_NonStrucAcc_BldgLvl + TotCompLoss_NonStrucAcc_FlrLvl;

    % Calculate total component loss (before considering demolition)
    TotLoss_IM = TotCompLoss_Struc + TotCompLoss_NonStrucDrift + TotCompLoss_NonStrucAcc;

    % =====================================================================
    % DEMOLITION LOSS ASSESSMENT
    % =====================================================================
    % Calculate expected loss due to building demolition based on
    % residual story drift ratios exceeding acceptable limits

    % Define residual drift range for PDF calculation
    x_RSDR_pdf = (0.0002:0.0002:0.06)';

    % Calculate maximum residual story drift for each ground motion
    Max_RSDR = zeros(1, length(RSDR));
    for gm = 1:length(RSDR)
        Max_RSDR(1, gm) = max(abs(RSDR(:, gm)));    % Maximum across all stories
    end

    % Demolition fragility curve parameters from Ramirez and Miranda (2012)
    % Based on engineering judgment and post-earthquake observations
    theta_hat_Demolish = 0.015;     % Median residual drift capacity (1.5%)
    beta_hat_Demolish = 0.3;        % Dispersion (log-standard deviation)

    % Calculate probability of demolition given residual drift
    % P(Demolition | RSDR) - uses normal CDF in log space
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

    % Transpose for consistency with other arrays
    PD_IM_NC = PD_IM_NC';

    % =====================================================================
    % REPLACEMENT COST CALCULATION
    % =====================================================================
    % Calculate building replacement cost based on floor area and unit cost

    Bay_width = Geometry(1, 6);         % Bay width [ft]
    Num_bays = 5;                       % Number of bays in each direction
    Floor_area = (Bay_width * Num_bays) * (Bay_width * Num_bays); % Floor area [ft²]
    Per_sqft_Cost = 250;                % Unit cost [$/ft²] from Hwang and Lignos (2017)
    Tot_Replac_Cost = Num_Storey * Floor_area * Per_sqft_Cost;   % Total replacement cost [$]

    % Demolition cost (10% of replacement cost for debris removal, etc.)
    DemoLoss = 0.1 * Tot_Replac_Cost;
    Tot_Loss_Demo = Tot_Replac_Cost + DemoLoss;  % Total loss if demolished

    % =====================================================================
    % EXPECTED LOSS CALCULATION WITH DEMOLITION CONSIDERATION
    % =====================================================================
    % Apply the demolition probability to calculate final expected losses
    % Loss = Repair_Loss * P(No_Demolition) + Demolition_Loss * P(Demolition)

    % Structural repair loss (conditional on no demolition)
    Loss_Struc(1, :) = TotCompLoss_Struc(1, :) .* (1. - PD_IM_NC(1, :));
    Norm_Loss_Struc(1, :) = Loss_Struc(1, :) ./ Tot_Replac_Cost;

    % Non-structural drift-sensitive repair loss (conditional on no demolition)
    Loss_NonStruc_Drift(1, :) = TotCompLoss_NonStrucDrift(1, :) .* (1. - PD_IM_NC(1, :));
    Norm_Loss_NonStruc_Drift(1, :) = Loss_NonStruc_Drift(1, :) ./ Tot_Replac_Cost;

    % Non-structural acceleration-sensitive repair loss (conditional on no demolition)
    Loss_NonStruc_Acc(1, :) = TotCompLoss_NonStrucAcc(1, :) .* (1. - PD_IM_NC(1, :));
    Norm_Loss_NonStruc_Acc(1, :) = Loss_NonStruc_Acc(1, :) ./ Tot_Replac_Cost;

    % Demolition loss (conditional on demolition decision)
    Loss_Demo(1, :) = Tot_Loss_Demo * (PD_IM_NC(1, :));
    Norm_Loss_Demo(1, :) = Loss_Demo(1, :) ./ Tot_Replac_Cost;

    % Total expected loss (repair + demolition)
    Loss_Tot(1, :) = TotLoss_IM(1, :) .* (1. - PD_IM_NC(1, :)) + Tot_Loss_Demo * (PD_IM_NC(1, :));
    Norm_Tot_Loss(1, :) = Loss_Tot(1, :) ./ Tot_Replac_Cost;

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
% excelFileName = 'FEMAP58_EAL_NormEAL_Latest-07-08-24.xlsx';
%
% % Specify the sheet name
% sheetName = 'Sheet1';
%
% % Write the data to Excel
% xlswrite(excelFileName, All_Output, sheetName);
% disp('Data written to Excel successfully.');


% Display the elapsed time in multiple units for performance monitoring
disp('=================================================================');
disp('                    FEMA P-58 ANALYSIS COMPLETE                 ');
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
% END OF FEMA P-58 SEISMIC LOSS ASSESSMENT
% =========================================================================
