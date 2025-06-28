%% Building Loss Assessment Using Story Loss Functions (SLFs)
% ========================================================================
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
% This MATLAB script performs comprehensive building loss assessment using
% pre-generated Story Loss Functions (SLFs). It computes Expected Annual 
% Loss (EAL) for structural and non-structural components considering
% demolition probability based on residual story drift.

% STORY LOSS FUNCTION (SLF) METHODOLOGY - THEORETICAL BACKGROUND
% ===============================================================
%
% Story Loss Functions (SLFs) are continuous functions that directly relate
% Engineering Demand Parameters (EDPs) to expected monetary losses at the
% story level, incorporating the effects of multiple components simultaneously.
%
% Reference:
% ----------
% Shahnazaryan, D., O'Reilly, G.J., and Monteiro, R. (2021). "Story loss 
% functions for seismic design and assessment: Development of tools and 
% application." Earthquake Spectra, 37(4), 2813-2839.
%
% INPUTS REQUIRED:
% ---------------
% 1. Ground Motion Intensity Measures: 'GuanDataBase-IMs.csv'
% 2. Structural Response Data: PSDR, PFA, RSDR files per building
% 3. Pre-generated Story Loss Functions: Structural, NSD, NSA components
% 4. Seismic Hazard Data: 'USGSHazard_data.csv'
%
% OUTPUTS GENERATED:
% -----------------
% - EAL values for each component type (structural, NSD, NSA, demolition)
% - Normalized EAL values (with respect to replacement cost)
% - Total expected annual loss per building
%
% HOW STORY LOSS FUNCTIONS ARE USED:
% ----------------------------------
% 1. SLF Selection: Based on building geometry (9 structural, 3 NSD, 15 NSA cases)
% 2. Demand-Loss Relationship: Fit regression, compute median demand, generate PDF
% 3. Loss Computation: Interpolate SLF, multiply by PDF, integrate over EDP range
% 4. EAL Integration: Weight by hazard curve derivative, integrate over IM range
%
% ========================================================================

clc; 
clearvars; close all;
clear;
tic;

%% Input GM IMs
% Load ground motion intensity measure data for all buildings
IM_FilePath = 'GuanDataBase-IMs.csv';

% Read the specified column from the file
Data = readmatrix(IM_FilePath);
IM_GM = Data(:, 3:242);   %%IM values for all the 240 ground motions buildingwise
ID = Data(:, 1);
IM_GM_trans = IM_GM';
Target_time_period = Data(:,2); 

% Define spectral acceleration values for loss computation
SaCalc=[0.001,0.01,0.05,0.1:0.01:2]';
% SaCalc=[0.1:0.01:2]';

%% Buildingwise Component Loss Assessment

% Initializing output arrays for Expected Annual Loss (EAL) computation
dlambdaSa=zeros(length(SaCalc),length(ID));
EAL_Tot= zeros(length(ID),1);
EAL_Demo= zeros(length(ID),1);
EAL_Struc = zeros(length(ID),1);
EAL_NonStrucDrift= zeros(length(ID),1);
EAL_NonStrucAcc = zeros(length(ID),1);       

% Normalized EAL arrays (with respect to replacement cost)
Norm_EAL_Tot= zeros(length(ID),1);
Norm_EAL_Demo = zeros(length(ID),1);
Norm_EAL_Struc = zeros(length(ID),1);
Norm_EAL_NonStrucDrift= zeros(length(ID),1);
Norm_EAL_NonStrucAcc = zeros(length(ID),1);

%%
% Main loop through all buildings for loss assessment
for bd = 1:length(ID)
    %% Reading necessary files
    % Define file paths for building-specific structural response data
    filePSDR = strcat(['Guan_database\Structura' ...
        'lResponses\EDPsUnder240GMs\PeakStoryDrift'], '\', 'Building_',num2str(ID(bd),'%.1d'), '_PSDR.csv');
    filePFA = strcat(['Guan_database\Structural' ...
        'Responses\EDPsUnder240GMs\PeakFloorAcceleration'], '\', 'Building_',num2str(ID(bd),'%.1d'), '_PFA.csv');
    fileRSDR = strcat(['Guan_database\Structural' ...
        'Responses\EDPsUnder240GMs\ResidualStoryDrift'], '\', 'Building_',num2str(ID(bd),'%.1d'), '_RSDR.csv');
    fileGeometry = strcat('Guan_database\BuildingDesigns', '\', 'Building_',num2str(ID(bd),'%.1d'), '\', 'Geometry.csv');
    
    % Load structural response data
    PSDR = readmatrix(filePSDR);
    PFA = readmatrix(filePFA);
    
    PFA=PFA*0.0254*0.10197; %% Convert PFA to g units
    RSDR = readmatrix(fileRSDR);
    
    PSDR_tans = PSDR';
    Geometry=readmatrix(fileGeometry); %% Read geometry file of Guan database
    
    %%  Structural Component Loss Assessment

    % Load pre-generated Story Loss Functions for structural components
    fileSLF_Struct_ColBasePlate = strcat('StoryLossFunctions_20000_old\MyInventory_Structural\MyInventory_Structural_Columnbase','\', 'CombinedOutput_SLF_Struct_ColBasePlate.csv'); %% Read fragility and cost data
    SLF_Struct_ColBasePlate = readmatrix(fileSLF_Struct_ColBasePlate);

    fileSLF_Struct_Splices = strcat('StoryLossFunctions_20000_old\MyInventory_Structural\MyInventory_Structural_Splices','\', 'CombinedOutput_SLF_Struct_Splices.csv'); %% Read fragility and cost data
    SLF_Struct_Splices = readmatrix(fileSLF_Struct_Splices);

    fileSLF_Struct_4comps = strcat('StoryLossFunctions_20000_old\MyInventory_Structural\MyInventory_Structural_4comp','\', 'CombinedOutput_SLF_Struct_4comps.csv'); %% Read fragility and cost data
    SLF_Struct_4comps = readmatrix(fileSLF_Struct_4comps);

    %%%%% Defining SLFs of structural components for each floor
    %%% Assigning SLFs based on Num_LFRS_Bays and bay-width (9 cases)
    % Extract building geometry parameters from the database
    Num_LFRS_Bays=Geometry(1,2);
    Num_stories=Geometry(1,1);
    Bay_width=Geometry(1,6);
    IDR_values=SLF_Struct_4comps(:,1);
    % IDR_values = round(IDR_values, 3);
    SLF_Splices=SLF_Struct_Splices(:,2);

    % Select appropriate Story Loss Functions based on building configuration
    % 9 different cases based on number of LFRS bays (1,3,5) and bay width (20,30,40 ft)
    if (Num_LFRS_Bays==1)&&(Bay_width==20)
        SLF_4COMP=SLF_Struct_4comps(:,2);
        SLF_ColBasePlate=SLF_Struct_ColBasePlate(:,2);
    elseif (Num_LFRS_Bays==1)&&(Bay_width==30)
        SLF_4COMP=SLF_Struct_4comps(:,3);
        SLF_ColBasePlate=SLF_Struct_ColBasePlate(:,2);
    elseif (Num_LFRS_Bays==1)&&(Bay_width==40)
        SLF_4COMP=SLF_Struct_4comps(:,4);
        SLF_ColBasePlate=SLF_Struct_ColBasePlate(:,2);
    elseif (Num_LFRS_Bays==3)&&(Bay_width==20)
        SLF_4COMP=SLF_Struct_4comps(:,5);
        SLF_ColBasePlate=SLF_Struct_ColBasePlate(:,5);
    elseif (Num_LFRS_Bays==3)&&(Bay_width==30)
        SLF_4COMP=SLF_Struct_4comps(:,6);
        SLF_ColBasePlate=SLF_Struct_ColBasePlate(:,5);
    elseif (Num_LFRS_Bays==3)&&(Bay_width==40)
        SLF_4COMP=SLF_Struct_4comps(:,7);
        SLF_ColBasePlate=SLF_Struct_ColBasePlate(:,5);
    elseif (Num_LFRS_Bays==5)&&(Bay_width==20)
        SLF_4COMP=SLF_Struct_4comps(:,8);
        SLF_ColBasePlate=SLF_Struct_ColBasePlate(:,8);
    elseif (Num_LFRS_Bays==5)&&(Bay_width==30)
        SLF_4COMP=SLF_Struct_4comps(:,9);
        SLF_ColBasePlate=SLF_Struct_ColBasePlate(:,8);
    elseif (Num_LFRS_Bays==5)&&(Bay_width==40)
        SLF_4COMP=SLF_Struct_4comps(:,10);
        SLF_ColBasePlate=SLF_Struct_ColBasePlate(:,8);
    end

    % Initialize arrays for structural component loss computation
    Loss_EDPSplices=zeros(height(PSDR),1);
    Loss_EDP4comp=zeros(height(PSDR),1);
    AllCompLoss_Struc=zeros(3,1);
    TotCompLoss_Struc=zeros(1,length(SaCalc));
    NumCompSets_Struct=3;  % Three component sets: base plates, splices, and four components

    % Loop through each intensity measure level for loss computation
    for im =1:length(SaCalc)
        % Process each structural component set
        for comp = 1:1:NumCompSets_Struct

            %% Component 1: Column Base Plates (Ground Floor Only)
            if comp==1
                flr=1;  % Ground floor only for base plates

                %fitting to Sa/MIDR - Develop demand model relating IM to EDP
                md_PIDR=fitlm(log(IM_GM(bd,:)),log(PSDR(flr,:)));
                %getting coefficients of regression
                PIDR_DemandPar(1:2)=md_PIDR.Coefficients{1:2,1};
                %getting sigma of the fitted model
                PIDR_DemandPar(3)=md_PIDR.RMSE;
                SIGMA = PIDR_DemandPar(3);
                MD_pidr = PIDR_DemandPar(2)*log(SaCalc(im))+PIDR_DemandPar(1);

                % Generate probability density function for IDR given IM
                PDF_PIDR_IM  = lognpdf(IDR_values,MD_pidr,SIGMA);

                Frag=PDF_PIDR_IM;

                MD_PIDR=exp(MD_pidr);  % Convert from log space


                %%%%%%% Enhanced interpolation approach for SLF evaluation
                if MD_PIDR < 0
                    Loss_EDP_Baseplate = 0;
                    % disp('MD_PIDR less than zero');
                else
                    if MD_PIDR <= min(IDR_values) || MD_PIDR >= max(IDR_values)
                        error('MD_PIDR is out of the range of IDR_values');
                    end

                    % Check if MD_PIDR is exactly equal to any value in IDR_values
                    found = false;
                    for i = 1:length(IDR_values)
                        if MD_PIDR == IDR_values(i)
                            Loss_EDP_Baseplate = SLF_ColBasePlate(i) .* Frag;
                            found = true;
                            break; % Exit the loop once the match is found
                        end
                    end

                    % If MD_PIDR is not exactly equal to any value, perform linear interpolation
                    if ~found
                        for i = 1:length(IDR_values) - 1
                            if MD_PIDR > IDR_values(i) && MD_PIDR < IDR_values(i + 1)
                                % Perform linear interpolation for SLF values
                                x1 = IDR_values(i);
                                x2 = IDR_values(i + 1);
                                y1_SLF_ColBasePlate = SLF_ColBasePlate(i);
                                y2_SLF_ColBasePlate = SLF_ColBasePlate(i + 1);
                                y1_Frag = Frag(i);
                                y2_Frag = Frag(i + 1);

                                % Interpolate SLF_ColBasePlate
                                SLF_ColBasePlate_interp = y1_SLF_ColBasePlate + (MD_PIDR - x1) * (y2_SLF_ColBasePlate - y1_SLF_ColBasePlate) / (x2 - x1);

                                % % Interpolate Frag (if needed)
                                % Frag_interp = y1_Frag + (MD_PIDR - x1) * (y2_Frag - y1_Frag) / (x2 - x1);

                                % Compute Loss_EDP_Baseplate
                                Loss_EDP_Baseplate = SLF_ColBasePlate_interp .* Frag;

                                break; % Exit the loop once the match is found
                            end
                        end
                    end
                end

                %%%%%%%
                % Integrate loss over the EDP range using trapezoidal rule
                CompLoss_EDP=trapz(IDR_values, Loss_EDP_Baseplate);
            end

            %% Component 2: Column Splices (Every Other Floor, Excluding Single-Story)
            if comp==2
                if Num_stories==1
                    CompLoss_EDP=0;  % No splices in single-story buildings
                else
                    for flr=1:2:height(PSDR)  % Every other floor
                        %fitting to Sa/MIDR - Develop demand model
                        md_PIDR=fitlm(log(IM_GM(bd,:)),log(PSDR(flr,:)));
                        %getting coefficients of regression
                        PIDR_DemandPar(1:2)=md_PIDR.Coefficients{1:2,1};
                        %getting sigma of the fitted model
                        PIDR_DemandPar(3)=md_PIDR.RMSE;
                        SIGMA = PIDR_DemandPar(3);
                        MD_pidr = PIDR_DemandPar(2)*log(SaCalc(im))+PIDR_DemandPar(1);

                        % Generate PDF for IDR given IM
                        PDF_PIDR_IM  = lognpdf(IDR_values,MD_pidr,SIGMA);

                        Frag=PDF_PIDR_IM;
                        MD_PIDR=exp(MD_pidr);

                        % SLF interpolation for splice components
                        if MD_PIDR < 0
                            Loss_EDP = 0;
                            disp('MD_PIDR less than zero');
                        else
                            if MD_PIDR <= min(IDR_values) || MD_PIDR >= max(IDR_values)
                                error('MD_PIDR is out of the range of IDR_values');
                            end

                            % Check if MD_PIDR is exactly equal to any value in IDR_values
                            found = false;
                            for i = 1:length(IDR_values)
                                if MD_PIDR == IDR_values(i)
                                    Loss_EDP = SLF_Splices(i) .* Frag;
                                    found = true;
                                    break; % Exit the loop once the match is found
                                end
                            end

                            % If MD_PIDR is not exactly equal to any value, perform linear interpolation
                            if ~found
                                for i = 1:length(IDR_values) - 1
                                    if MD_PIDR > IDR_values(i) && MD_PIDR < IDR_values(i + 1)
                                        % Perform linear interpolation
                                        x1 = IDR_values(i);
                                        x2 = IDR_values(i + 1);
                                        y1_SLF_Splices = SLF_Splices(i);
                                        y2_SLF_Splices = SLF_Splices(i + 1);
                                        y1_Frag = Frag(i);
                                        y2_Frag = Frag(i + 1);

                                        % Interpolate SLF_Splices
                                        SLF_Splices_interp = y1_SLF_Splices + (MD_PIDR - x1) * (y2_SLF_Splices - y1_SLF_Splices) / (x2 - x1);

                                        % % Interpolate Frag
                                        % Frag_interp = y1_Frag + (MD_PIDR - x1) * (y2_Frag - y1_Frag) / (x2 - x1);

                                        % Compute Loss_EDP
                                        Loss_EDP = SLF_Splices_interp .* Frag;

                                        break;
                                    end
                                end
                            end
                        end

                        Loss_EDPSplices(flr,1) = trapz(IDR_values, Loss_EDP); %%%floor level loss
                    end
                    Loss_EDP_Splices=sum(Loss_EDPSplices(:,1));
                    CompLoss_EDP=Loss_EDP_Splices;
                end
            end

            %% Component 3: Remaining Four Components (All Floors)
            % Includes shear tabs, RBS connections, and slab components
            if comp==3
                for flr=1:1:height(PSDR)  % All floors
                    %fitting to Sa/MIDR - Develop demand model
                    md_PIDR=fitlm(log(IM_GM(bd,:)),log(PSDR(flr,:)));
                    %getting coefficients of regression
                    PIDR_DemandPar(1:2)=md_PIDR.Coefficients{1:2,1};
                    %getting sigma of the fitted model
                    PIDR_DemandPar(3)=md_PIDR.RMSE;
                    SIGMA = PIDR_DemandPar(3);
                    MD_pidr = PIDR_DemandPar(2)*log(SaCalc(im))+PIDR_DemandPar(1);

                    % Generate PDF for IDR given IM
                    PDF_PIDR_IM  = lognpdf(IDR_values,MD_pidr,SIGMA);

                    Frag=PDF_PIDR_IM;

                    MD_PIDR=exp(MD_pidr);

                    % Enhanced interpolation approach for four-component SLF
                    if MD_PIDR < 0
                        Loss_EDP = 0;
                        disp('MD_PIDR less than zero');
                    else
                        if MD_PIDR <= min(IDR_values) || MD_PIDR >= max(IDR_values)
                            error('MD_PIDR is out of the range of IDR_values');
                        end

                        % Check if MD_PIDR is exactly equal to any value in IDR_values
                        found = false;
                        for i = 1:length(IDR_values)
                            if MD_PIDR == IDR_values(i)
                                Loss_EDP = SLF_4COMP(i) .* Frag;
                                found = true;
                                break; % Exit the loop once the match is found
                            end
                        end

                        % If MD_PIDR is not exactly equal to any value, perform linear interpolation
                        if ~found
                            for i = 1:length(IDR_values) - 1
                                if MD_PIDR > IDR_values(i) && MD_PIDR < IDR_values(i + 1)
                                    % Perform linear interpolation
                                    x1 = IDR_values(i);
                                    x2 = IDR_values(i + 1);
                                    y1_SLF_4COMP = SLF_4COMP(i);
                                    y2_SLF_4COMP = SLF_4COMP(i + 1);
                                    y1_Frag = Frag(i);
                                    y2_Frag = Frag(i + 1);

                                    % Interpolate SLF_4COMP
                                    SLF_4COMP_interp = y1_SLF_4COMP + (MD_PIDR - x1) * (y2_SLF_4COMP - y1_SLF_4COMP) / (x2 - x1);

                                    % % Interpolate Frag
                                    % Frag_interp = y1_Frag + (MD_PIDR - x1) * (y2_Frag - y1_Frag) / (x2 - x1);

                                    % Compute Loss_EDP
                                    Loss_EDP = SLF_4COMP_interp .* Frag;

                                    break; % Exit the loop once the match is found
                                end
                            end
                        end
                    end

                    Loss_EDP4comp(flr,1) = trapz(IDR_values, Loss_EDP); %%%floor level loss
                end
                Loss_EDP_4comps=sum(Loss_EDP4comp(:,1));
                CompLoss_EDP=Loss_EDP_4comps;
            end

            AllCompLoss_Struc(comp,1)=CompLoss_EDP;
        end
        % Sum losses from all structural component sets
        TotCompLoss_Struc(1,im)=sum(AllCompLoss_Struc(:,1));
    end
    
    %% Non-Structural Drift-Sensitive Loss Assessment

    % Load SLFs for non-structural drift-sensitive components
    fileSLF_NSD = strcat('StoryLossFunctions_20000_old\MyInventory_NSD','\', 'CombinedOutput_SLF_NSD.csv');
    SLF_NSD_ = readmatrix(fileSLF_NSD);

    %%% Assigning SLFs based on bay-width (3 cases)
    Bay_width=Geometry(1,6);

    IDR_values=SLF_NSD_(:,1);
    % IDR_values = round(IDR_values, 3);

    % Select appropriate NSD SLF based on bay width
    if Bay_width==20
        SLF_NSD=SLF_NSD_(:,2);
    elseif Bay_width==30
        SLF_NSD=SLF_NSD_(:,3);
    elseif Bay_width==40
        SLF_NSD=SLF_NSD_(:,4);
    end

    Loss_EDP_NSD=zeros(height(PSDR),1);
    TotCompLoss_NonStrucDrift=zeros(1,length(SaCalc));

    % Loop through each IM level and floor for NSD components
    for im =1:length(SaCalc)

        for flr=1:1:height(PSDR)

            %%%%%%%%%%%%%%%%%%%% Computing Loss %%%%%%%%%%%%%%%%%%%%%%%

            % Develop demand model for NSD components
            md_PIDR=fitlm(log(IM_GM(bd,:)),log(PSDR(flr,:)));
            %getting coefficients of regression
            PIDR_DemandPar(1:2)=md_PIDR.Coefficients{1:2,1};
            %getting sigma of the fitted model
            PIDR_DemandPar(3)=md_PIDR.RMSE;
            SIGMA = PIDR_DemandPar(3);
            MD_pidr = PIDR_DemandPar(2)*log(SaCalc(im))+PIDR_DemandPar(1);

            % Generate PDF for IDR given IM
            PDF_PIDR_IM  = lognpdf(IDR_values,MD_pidr,SIGMA);

            Frag=PDF_PIDR_IM;

            MD_PIDR=exp(MD_pidr);
            % Round MD_PIDR to three decimal places
            % MD_PIDR = round(MD_PIDR, 3);

            % SLF interpolation for NSD components
            if MD_PIDR < 0
                Loss_EDP = 0;
                disp('MD_PIDR less than zero');
            else
                if MD_PIDR <= min(IDR_values) || MD_PIDR >= max(IDR_values)
                    error('MD_PIDR is out of the range of IDR_values');
                end

                % Check if MD_PIDR is exactly equal to any value in IDR_values
                found = false;
                for i = 1:length(IDR_values)
                    if MD_PIDR == IDR_values(i)
                        Loss_EDP = SLF_NSD(i) .* Frag;
                        found = true;
                        break; % Exit the loop once the match is found
                    end
                end

                % If MD_PIDR is not exactly equal to any value, perform linear interpolation
                if ~found
                    for i = 1:length(IDR_values) - 1
                        if MD_PIDR > IDR_values(i) && MD_PIDR < IDR_values(i + 1)
                            % Perform linear interpolation
                            x1 = IDR_values(i);
                            x2 = IDR_values(i + 1);
                            y1_SLF_NSD = SLF_NSD(i);
                            y2_SLF_NSD = SLF_NSD(i + 1);
                            y1_Frag = Frag(i);
                            y2_Frag = Frag(i + 1);

                            % Interpolate SLF_NSD
                            SLF_NSD_interp = y1_SLF_NSD + (MD_PIDR - x1) * (y2_SLF_NSD - y1_SLF_NSD) / (x2 - x1);

                            % % Interpolate Frag
                            % Frag_interp = y1_Frag + (MD_PIDR - x1) * (y2_Frag - y1_Frag) / (x2 - x1);

                            % Compute Loss_EDP
                            Loss_EDP = SLF_NSD_interp .* Frag;

                            break; % Exit the loop once the match is found
                        end
                    end
                end
            end

            Loss_EDP_NSD(flr,1) = trapz(IDR_values, Loss_EDP); %%%floor level loss
        end

        TotCompLoss_NonStrucDrift(1,im)=sum(Loss_EDP_NSD(:,1));
    end

    %% Non-Structural Acceleration-Sensitive Loss Assessment

    %%%%%%%%%%%%%%%%%%%%%%% Floor Level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Load SLFs for floor-level non-structural acceleration-sensitive components
    fileSLF_NSA_FLRLVL = strcat('StoryLossFunctions_20000_old\MyInventory_NSA\MyInventory_NSA_FloorLevelComponents','\', 'CombinedOutput_SLF_NSA_Flrlvl.csv');
    SLF_NSA_FLRLVL = readmatrix(fileSLF_NSA_FLRLVL);

    %%% Assigning SLFs based on bay-width (3 cases)
    Bay_width=Geometry(1,6);

    PFA_values=SLF_NSA_FLRLVL(:,1);
    % PFA_values = round(PFA_values, 3);

    % Select appropriate floor-level NSA SLF based on bay width
    if Bay_width==20
        SLF_NSA_flrlvl=SLF_NSA_FLRLVL(:,2);
    elseif Bay_width==30
        SLF_NSA_flrlvl=SLF_NSA_FLRLVL(:,3);
    elseif Bay_width==40
        SLF_NSA_flrlvl=SLF_NSA_FLRLVL(:,4);
    end

    Loss_EDP_NSA_FlrLvl=zeros(height(PSDR),1);
    TotCompLoss_NonStrucAcc_FlrLvl=zeros(1,length(SaCalc));

    % Process floor-level NSA components
    for im =1:length(SaCalc)

        for flr=1:1:height(PFA)-1

            % Develop demand model for Peak Floor Acceleration (PFA)
            md_PFA=fitlm(log(IM_GM(bd,:)),log(PFA(flr,:)));
            %getting coefficients of regression
            PFA_DemandPar(1:2)=md_PFA.Coefficients{1:2,1};
            %getting sigma of the fitted model
            PFA_DemandPar(3)=md_PFA.RMSE;
            SIGMA = PFA_DemandPar(3);
            MD_pfa = PFA_DemandPar(2)*log(SaCalc(im))+PFA_DemandPar(1);

            % Generate PDF for PFA given IM
            PDF_PFA_IM  = lognpdf(PFA_values,MD_pfa,SIGMA);
            Frag=PDF_PFA_IM;

            MD_PFA=exp(MD_pfa);
            % SLF interpolation for floor-level NSA components
            if MD_PFA < 0
                Loss_EDP = 0;
                disp('MD_PFA less than zero');
            else
                if MD_PFA <= min(PFA_values) || MD_PFA >= max(PFA_values)
                    error('MD_PFA is out of the range of PFA_values');
                end

                % Check if MD_PFA is exactly equal to any value in PFA_values
                found = false;
                for i = 1:length(PFA_values)
                    if MD_PFA == PFA_values(i)
                        Loss_EDP = SLF_NSA_flrlvl(i) .* Frag;
                        found = true;
                        break; % Exit the loop once the match is found
                    end
                end

                % If MD_PFA is not exactly equal to any value, perform linear interpolation
                if ~found
                    for i = 1:length(PFA_values) - 1
                        if MD_PFA > PFA_values(i) && MD_PFA < PFA_values(i + 1)
                            % Perform linear interpolation
                            x1 = PFA_values(i);
                            x2 = PFA_values(i + 1);
                            y1_SLF_NSA_flrlvl = SLF_NSA_flrlvl(i);
                            y2_SLF_NSA_flrlvl = SLF_NSA_flrlvl(i + 1);
                            y1_Frag = Frag(i);
                            y2_Frag = Frag(i + 1);

                            % Interpolate SLF_NSA_flrlvl
                            SLF_NSA_flrlvl_interp = y1_SLF_NSA_flrlvl + (MD_PFA - x1) * (y2_SLF_NSA_flrlvl - y1_SLF_NSA_flrlvl) / (x2 - x1);

                            % % Interpolate Frag
                            % Frag_interp = y1_Frag + (MD_PFA - x1) * (y2_Frag - y1_Frag) / (x2 - x1);

                            % Compute Loss_EDP
                            Loss_EDP = SLF_NSA_flrlvl_interp .* Frag;

                            break; % Exit the loop once the match is found
                        end
                    end
                end
            end

            Loss_EDP_NSA_FlrLvl(flr,1) = trapz(PFA_values, Loss_EDP); %%%floor level loss
        end
        TotCompLoss_NonStrucAcc_FlrLvl(1,im)=sum(Loss_EDP_NSA_FlrLvl(:,1));
    end
   
    %%%%%%%%%%%%%%%%%%%%%% Building Level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Load SLFs for building-level non-structural acceleration-sensitive components
    fileSLF_NSA_BldgLvl = strcat('StoryLossFunctions_20000_old\MyInventory_NSA\MyInventory_NSA_BuildingLevelComponents','\', 'CombinedOutput_SLF_NSA_Bldglvl.csv');
    SLF_NSA_BldgLvl_ = readmatrix(fileSLF_NSA_BldgLvl);

    %%% Assigning component quantities based on number of stories
    %%% and bay-width (15 cases)
    Bay_width=Geometry(1,6);
    Num_Storey=Geometry(1,1);

    PFA_values=SLF_NSA_BldgLvl_(:,1);
    % PFA_values = round(PFA_values, 3);

    % Select appropriate building-level NSA SLF based on building configuration
    % 15 different cases based on number of stories (1,5,9,14,19) and bay width (20,30,40 ft)
    if (Num_Storey==1)&&(Bay_width==20)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,2);
    elseif (Num_Storey==1)&&(Bay_width==30)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,3);
    elseif (Num_Storey==1)&&(Bay_width==40)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,4);

    elseif (Num_Storey==5)&&(Bay_width==20)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,5);
    elseif (Num_Storey==5)&&(Bay_width==30)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,6);
    elseif (Num_Storey==5)&&(Bay_width==40)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,7);

    elseif (Num_Storey==9)&&(Bay_width==20)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,8);
    elseif (Num_Storey==9)&&(Bay_width==30)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,9);
    elseif (Num_Storey==9)&&(Bay_width==40)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,10);

    elseif (Num_Storey==14)&&(Bay_width==20)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,11);
    elseif (Num_Storey==14)&&(Bay_width==30)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,12);
    elseif (Num_Storey==14)&&(Bay_width==40)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,13);

    elseif (Num_Storey==19)&&(Bay_width==20)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,14);
    elseif (Num_Storey==19)&&(Bay_width==30)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,15);
    elseif (Num_Storey==19)&&(Bay_width==40)
        SLF_NSA_BldgLvl=SLF_NSA_BldgLvl_(:,16);
    end

    Loss_EDP_NSA_BldgLvl=zeros(height(PSDR),1);
    TotCompLoss_NonStrucAcc_BldgLvl =zeros(1,length(SaCalc));

    % Process building-level NSA components
    for im =1:length(SaCalc)
        flr=1;  % Use ground floor PFA for building-level components

        % Develop demand model for building-level NSA components
        md_PFA=fitlm(log(IM_GM(bd,:)),log(PFA(flr,:)));
        %getting coefficients of regression
        PFA_DemandPar(1:2)=md_PFA.Coefficients{1:2,1};
        %getting sigma of the fitted model
        PFA_DemandPar(3)=md_PFA.RMSE;
        SIGMA = PFA_DemandPar(3);
        MD_pfa = PFA_DemandPar(2)*log(SaCalc(im))+PFA_DemandPar(1);

        % Generate PDF for PFA given IM
        PDF_PFA_IM  = lognpdf(PFA_values,MD_pfa,SIGMA);
        Frag=PDF_PFA_IM;

        MD_PFA=exp(MD_pfa);
        % Round MD_PIDR to three decimal places
        % MD_PFA = round(MD_PFA, 3);

        % SLF interpolation for building-level NSA components
        if MD_PFA < 0
            Loss_EDP = 0;
            disp('MD_PFA less than zero');
        else
            if MD_PFA <= min(PFA_values) || MD_PFA >= max(PFA_values)
                error('MD_PFA is out of the range of PFA_values');
            end

            % Check if MD_PFA is exactly equal to any value in PFA_values
            found = false;
            for i = 1:length(PFA_values)
                if MD_PFA == PFA_values(i)
                    Loss_EDP = SLF_NSA_BldgLvl(i) .* Frag;
                    found = true;
                    break; % Exit the loop once the match is found
                end
            end

            % If MD_PFA is not exactly equal to any value, perform linear interpolation
            if ~found
                for i = 1:length(PFA_values) - 1
                    if MD_PFA > PFA_values(i) && MD_PFA < PFA_values(i + 1)
                        % Perform linear interpolation
                        x1 = PFA_values(i);
                        x2 = PFA_values(i + 1);
                        y1_SLF_NSA_BldgLvl = SLF_NSA_BldgLvl(i);
                        y2_SLF_NSA_BldgLvl = SLF_NSA_BldgLvl(i + 1);
                        y1_Frag = Frag(i);
                        y2_Frag = Frag(i + 1);

                        % Interpolate SLF_NSA_BldgLvl
                        SLF_NSA_BldgLvl_interp = y1_SLF_NSA_BldgLvl + (MD_PFA - x1) * (y2_SLF_NSA_BldgLvl - y1_SLF_NSA_BldgLvl) / (x2 - x1);

                        % % Interpolate Frag
                        % Frag_interp = y1_Frag + (MD_PFA - x1) * (y2_Frag - y1_Frag) / (x2 - x1);

                        % Compute Loss_EDP
                        Loss_EDP = SLF_NSA_BldgLvl_interp .* Frag;

                        break; % Exit the loop once the match is found
                    end
                end
            end
        end
        Loss_EDP_NSA_BldgLvl(flr,1) = trapz(PFA_values, Loss_EDP); %%%floor level loss
        TotCompLoss_NonStrucAcc_BldgLvl(1,im)=sum(Loss_EDP_NSA_BldgLvl(:,1));
    end

    % Combine floor-level and building-level NSA losses
    TotCompLoss_NonStrucAcc=TotCompLoss_NonStrucAcc_BldgLvl+TotCompLoss_NonStrucAcc_FlrLvl;

    % Total loss from all components (before considering demolition)
    TotLoss_IM=TotCompLoss_Struc+TotCompLoss_NonStrucDrift+TotCompLoss_NonStrucAcc;

    %% Expected loss due to demolition (Demolition Loss)

    %%% Compute probability that the building is being considered to be demolished
    % First, compute probability density function of maximum residual story drift ratio, fRSDR|IM
    x_RSDR_pdf = (0.0002:0.0002:0.06)';
    Max_RSDR=zeros(1,length(RSDR));

    % Find maximum RSDR across all floors for each ground motion
    for gm=1:1:length(RSDR)
        Max_RSDR(1,gm)=max(abs(RSDR(:,gm)));
    end

    % Assumed fragility curve for decision of demolition of building based
    % on Ramirez and Miranda (2012)...
    theta_hat_Demolish = 0.015; beta_hat_Demolish = 0.3; % From Ramirez and Miranda (2012)
    PD_RSDR = normcdf((log(x_RSDR_pdf/theta_hat_Demolish))/beta_hat_Demolish);
    % CDF_RSDR = logncdf(x_RSDR_pdf, log(theta_hat_Demolish), beta_hat_Demolish);
    % PD_RSDR = [0; diff(CDF_RSDR)];
    % PD_RSDR = lognpdf(x_RSDR_pdf, log(theta_hat_Demolish), beta_hat_Demolish);


    PD_IM_NC = zeros(length(SaCalc),1);

    % Compute demolition probability for each IM level
    for im =1:length(SaCalc)

        % Develop demand model for maximum RSDR
        md_RSDR=fitlm(log(IM_GM(bd,:)),log(Max_RSDR(1,:)));
        %getting coefficients of regression
        RSDR_DemandPar(1:2)=md_RSDR.Coefficients{1:2,1};
        %getting sigma of the fitted model
        RSDR_DemandPar(3)=md_RSDR.RMSE;
        SIGMA_RSDR = RSDR_DemandPar(3);
        MD_RSDR = RSDR_DemandPar(2)*log(SaCalc(im))+RSDR_DemandPar(1);

        % PDF of RSDR at given IM level
        dP_RIDR_IM  = lognpdf(x_RSDR_pdf,MD_RSDR,SIGMA_RSDR);    % PDF of RSDR at IM=im

        % Convolution to get demolition probability
        PD_RSDR_x_dP_RSDR_IM = PD_RSDR(:,1) .* dP_RIDR_IM(:,1);

        PD_IM_NC(im,1) = trapz(x_RSDR_pdf, PD_RSDR_x_dP_RSDR_IM); % Probability that building is being considered to be demolished

    end
    PD_IM_NC = PD_IM_NC';

    %%%%%%%%%%%%%%%%%%%%%%%% Expected loss due to demolition and Replacement Cost %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate replacement cost and demolition loss
    Bay_width=Geometry(1,6);
    Num_bays=5;  % Assumed number of bays in each direction
    Floor_area= (Bay_width*Num_bays)*(Bay_width*Num_bays);  % Total floor area
    Per_sqft_Cost=250; %% 250 dollars per square foot adopted from Hwang and Lignos (2017)
    Tot_Replac_Cost=Num_Storey*Floor_area*Per_sqft_Cost;  % Total replacement cost
    DemoLoss=0.1*Tot_Replac_Cost; %% 10% of replacement cost for demolition activities
    Tot_Loss_Demo=Tot_Replac_Cost+DemoLoss;  % Total loss if building is demolished

    % ------------------------------------------------------------------------------------------
    % MULTIPLY PROBABILITY OF DEMOLISH
    % ------------------------------------------------------------------------------------------

    % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
    % COMPUTE: "Structural repair loss"  - - - - - - - - - - - - -
    % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
    % Structural repair loss (reduced by demolition probability)
    Loss_Struc(1,:) = TotCompLoss_Struc(1,:) .* (1.-PD_IM_NC(1,:));
    Norm_Loss_Struc(1,:)=Loss_Struc(1,:)./Tot_Replac_Cost;

    % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
    % COMPUTE: "Non-structural repair loss (drift)"  - - - - - - -
    % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
    % Non-structural drift-sensitive repair loss (reduced by demolition probability)
    Loss_NonStruc_Drift(1,:) = TotCompLoss_NonStrucDrift(1,:) .* (1.-PD_IM_NC(1,:));
    Norm_Loss_NonStruc_Drift(1,:)=Loss_NonStruc_Drift(1,:)./Tot_Replac_Cost;

    % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
    % COMPUTE: "Non-structural repair loss (acc)"  - - - - - - - -
    % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
    % Non-structural acceleration-sensitive repair loss (reduced by demolition probability)
    Loss_NonStruc_Acc(1,:) = TotCompLoss_NonStrucAcc(1,:)  .* (1.-PD_IM_NC(1,:));
    Norm_Loss_NonStruc_Acc(1,:)=Loss_NonStruc_Acc(1,:)./Tot_Replac_Cost;

    % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
    % COMPUTE: "Demolition loss"   - - - - - - - - - - - - - - - - -
    % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
    % Demolition loss (weighted by demolition probability)
    Loss_Demo(1,:) = Tot_Loss_Demo  * (PD_IM_NC(1,:));
    Norm_Loss_Demo(1,:)=Loss_Demo(1,:)./Tot_Replac_Cost;

    % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
    % COMPUTE: "Total loss"   - - - - - - - - - - - - - - - - -
    % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
    % Total loss (repair if not demolished + total replacement if demolished)
    Loss_Tot(1,:) = TotLoss_IM(1,:) .* (1.-PD_IM_NC(1,:))+ Tot_Loss_Demo  * (PD_IM_NC(1,:));
    Norm_Tot_Loss(1,:)=Loss_Tot(1,:)./Tot_Replac_Cost;

    %% Hazard Input
    % Load and process seismic hazard data
    Hazard_file = 'USGSHazard_data.csv';
    HazardData = readmatrix(Hazard_file);
    AllLamda_IM = HazardData(:,2:12);  % Hazard curves for different periods
    IM_values = HazardData(:,1);       % IM values for hazard curves

    % Extract hazard data for different periods
    a = AllLamda_IM(:,2);   % 0.1s
    b = AllLamda_IM(:,3);   % 0.2s
    c = AllLamda_IM(:,4);   % 0.3s
    d = AllLamda_IM(:,5);   % 0.5s
    e = AllLamda_IM(:,6);   % 0.75s
    f = AllLamda_IM(:,7);   % 1.0s
    g = AllLamda_IM(:,8);   % 2.0s
    h = AllLamda_IM(:,9);   % 3.0s
    i = AllLamda_IM(:,10);  % 4.0s
    j = AllLamda_IM(:,11);  % 5.0s

    Lamda_IM_target = zeros(length(IM_values),length(Target_time_period));

    % Interpolate hazard curve based on building's target period
    if Target_time_period(bd,1)>=0.1 && Target_time_period(bd,1)<0.2
        for nn = 1:1: length(IM_values)
            Lamda_IM_target(nn,bd) = interp1([0.1, 0.2], [log(a(nn,1)), log(b(nn,1))], Target_time_period(bd,1), 'linear');
        end

    elseif Target_time_period(bd,1)>=0.2 && Target_time_period(bd,1)<0.3
        for nn = 1:1: length(IM_values)
            Lamda_IM_target(nn,bd) = interp1([0.2, 0.3], [log(b(nn,1)), log(c(nn,1))], Target_time_period(bd,1), 'linear');
        end

    elseif Target_time_period(bd,1)>=0.3 && Target_time_period(bd,1)<0.5
        for nn = 1:1: length(IM_values)
            Lamda_IM_target(nn,bd) = interp1([0.3, 0.5], [log(c(nn,1)), log(d(nn,1))], Target_time_period(bd,1), 'linear');
        end

    elseif Target_time_period(bd,1)>=0.5 && Target_time_period(bd,1)<0.75
        for nn = 1:1: length(IM_values)
            Lamda_IM_target(nn,bd) = interp1([0.5, 0.75], [log(d(nn,1)), log(e(nn,1))], Target_time_period(bd,1), 'linear');
        end

    elseif Target_time_period(bd,1)>=0.75 && Target_time_period(bd,1)<1
        for nn = 1:length(IM_values)
            Lamda_IM_target(nn,bd) = interp1([0.75, 1], [log(e(nn,1)), log(f(nn,1))], Target_time_period(bd,1), 'linear');
        end

    elseif Target_time_period(bd,1)>=1 && Target_time_period(bd,1)<2
        for nn = 1:1: length(IM_values)
            Lamda_IM_target(nn,bd) = interp1([1, 2], [log(f(nn,1)), log(g(nn,1))], Target_time_period(bd,1), 'linear');
        end

    elseif Target_time_period(bd,1)>=2 && Target_time_period(bd,1)<3
        for nn = 1:1: length(IM_values)
            Lamda_IM_target(nn,bd) = interp1([2, 3], [log(g(nn,1)), log(h(nn,1))], Target_time_period(bd,1), 'linear');
        end

    elseif Target_time_period(bd,1)>=3 && Target_time_period(bd,1)<4
        for nn = 1:1: length(IM_values)
            Lamda_IM_target(nn,bd) = interp1([3, 4], [log(h(nn,1)), log(i(nn,1))], Target_time_period(bd,1), 'linear');
        end

    elseif Target_time_period(bd,1)>=4 && Target_time_period(bd,1)<5
        for nn = 1:1: length(IM_values)
            Lamda_IM_target(nn,bd) = interp1([4, 5], [log(i(nn,1)), log(j(nn,1))], Target_time_period(bd,1), 'linear');
        end
    end

    % Reformat as vector with SA (SaCalc) and exceedance rate
    hazardCurveSa=[IM_values,exp(Lamda_IM_target(:,bd))];
    PDF_temp=zeros(length(SaCalc),1);

    % Compute derivative of hazard curve using numerical differentiation
    hh=0.0000001;  % Small increment for numerical differentiation
    for kk = 1:length(SaCalc)
        % Interpolating value of F(x+hh)
        fxp = interp1(hazardCurveSa(:,1),hazardCurveSa(:,2),SaCalc(kk)+hh,'spline');
        %Interpolating value of F(x-hh)
        fxn = interp1(hazardCurveSa(:,1),hazardCurveSa(:,2),SaCalc(kk)-hh,'spline');
        % Two-point formula to compute numerical differentiation
        PDF_temp(kk) = abs((fxp - fxn))/(2*hh);
    end
    dlambdaSa(:,bd) = smooth(PDF_temp);  % Smooth for numerical stability

    %%   EAL Computation

    % Compute convolution of loss functions with hazard curve derivative for EAL integration
    Conv_loss_Struc = Loss_Struc(1,:)'.* dlambdaSa(:,bd);
    Conv_loss_NonStrucDrift = Loss_NonStruc_Drift(1,:)'.* dlambdaSa(:,bd);
    Conv_loss_NonStrucAcc = Loss_NonStruc_Acc(1,:)'.* dlambdaSa(:,bd);
    Conv_loss_Tot= Loss_Tot(1,:)'.* dlambdaSa(:,bd);
    Conv_loss_Demo= Loss_Demo(1,:)'.* dlambdaSa(:,bd);

    % Compute normalized convolution (with respect to replacement cost)
    Conv_Norm_loss_Struc = Norm_Loss_Struc(1,:)'.* dlambdaSa(:,bd);
    Conv_Norm_loss_NonStrucDrift = Norm_Loss_NonStruc_Drift(1,:)'.* dlambdaSa(:,bd);
    Conv_Norm_loss_NonStrucAcc = Norm_Loss_NonStruc_Acc(1,:)'.* dlambdaSa(:,bd);
    Conv_Norm_loss_Tot= Norm_Tot_Loss(1,:)'.* dlambdaSa(:,bd);
    Conv_Norm_loss_Demo= Norm_Loss_Demo(1,:)'.* dlambdaSa(:,bd);

    % Integrate over the IM range to compute Expected Annual Loss (EAL)
    EAL_Tot(bd,1) = trapz(SaCalc,Conv_loss_Tot);
    EAL_Demo(bd,1) = trapz(SaCalc,Conv_loss_Demo);
    EAL_Struc(bd,1) = trapz(SaCalc,Conv_loss_Struc);
    EAL_NonStrucDrift(bd,1) = trapz(SaCalc,Conv_loss_NonStrucDrift);
    EAL_NonStrucAcc(bd,1) = trapz(SaCalc,Conv_loss_NonStrucAcc);

    %%% Normalized EAL (with respect to Replacement cost)

    Norm_EAL_Tot(bd,1) = trapz(SaCalc,Conv_Norm_loss_Tot);
    Norm_EAL_Demo(bd,1) = trapz(SaCalc,Conv_Norm_loss_Demo);
    Norm_EAL_Struc(bd,1) = trapz(SaCalc,Conv_Norm_loss_Struc);
    Norm_EAL_NonStrucDrift(bd,1) = trapz(SaCalc,Conv_Norm_loss_NonStrucDrift);
    Norm_EAL_NonStrucAcc(bd,1) = trapz(SaCalc,Conv_Norm_loss_NonStrucAcc);

end  % End of main building loop

%% Compile output data matrix
% Organize all results into a single output matrix for easy analysis and export
All_Output(:,1)=ID(:,1);                        % Building IDs
All_Output(:,2)=EAL_Struc(:,1);                 % Structural EAL ($)
All_Output(:,3)=EAL_NonStrucDrift(:,1);         % Non-structural drift-sensitive EAL ($)
All_Output(:,4)=EAL_NonStrucAcc(:,1);           % Non-structural acceleration-sensitive EAL ($)
All_Output(:,5)=EAL_Demo(:,1);                  % Demolition EAL ($)
All_Output(:,6)=EAL_Tot(:,1);                   % Total EAL ($)
All_Output(:,7)=Norm_EAL_Struc(:,1);            % Normalized structural EAL (fraction of replacement cost)
All_Output(:,8)=Norm_EAL_NonStrucDrift(:,1);    % Normalized non-structural drift EAL (fraction)
All_Output(:,9)=Norm_EAL_NonStrucAcc(:,1);      % Normalized non-structural acceleration EAL (fraction)
All_Output(:,10)=Norm_EAL_Demo(:,1);            % Normalized demolition EAL (fraction)
All_Output(:,11)=Norm_EAL_Tot(:,1);             % Normalized total EAL (fraction)

% Optional: Specify the Excel file name for output export
% Uncomment the following lines to export results to Excel
% excelFileName = 'SLF_EAL_NormEAL_Results.xlsx';
% sheetName = 'EAL_Results';
% xlswrite(excelFileName, All_Output, sheetName);
% disp('Data written to Excel successfully.');

% Stop the timer and display performance statistics
elapsedTime = toc;

% Display the elapsed time and completion message
disp(['Elapsed Time: ' num2str(elapsedTime) ' seconds']);
disp(['Elapsed Time: ' num2str(elapsedTime / 60) ' minutes']);
disp(['Elapsed Time: ' num2str(elapsedTime / (60*60)) ' hours']);

% Display completion message with summary
fprintf('\n=== SLF-BASED SEISMIC LOSS ASSESSMENT COMPLETED ===\n');
fprintf('Total buildings processed: %d\n', length(ID));
fprintf('Average processing time per building: %.3f seconds\n', elapsedTime/length(ID));
fprintf('Output matrix size: %d x %d\n', size(All_Output,1), size(All_Output,2));
fprintf('Results ready for analysis and export.\n');
fprintf('=====================================================\n');
