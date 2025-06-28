README for scripts and data associated with the paper "Impact of Modeling Decisions on Seismic Loss and Fragility Assessment of Steel Buildings"

This repository contains scripts and data associated with the research paper "Impact of Modeling Decisions on Seismic Loss and Fragility Assessment of Steel Buildings." The scripts and data can be used to reproduce analyses and plots associated with the companion journal manuscript. This repository should be cited as:

**Citation:**
When using this repository, please cite:

```
Baddipalli, S., Gentile, R., O'Reilly, G.J., Shahnazaryan, D., & Esteghamati, M.Z. (2025). 
Impact of Modeling Decisions on Seismic Loss and Fragility Assessment of Steel Buildings. 
Earthquake Spectra. [Under Review]
```

**Repository Citation:**
```
Baddipalli, S., Gentile, R., O'Reilly, G.J., Shahnazaryan, D., & Esteghamati, M.Z. (2025). 
Scripts and data for Impact of Modeling Decisions on Seismic Loss and Fragility Assessment 
of Steel Buildings. GitHub Repository. https://github.com/yourusername/repository-name
```

## Table of Contents
- [1. Overview and Motivation](#1-overview-and-motivation)
- [2. Methodology Overview](#2-methodology-overview)
  * [2.1 Methodology Flowchart](#21-methodology-flowchart)
  * [2.2 Database Description](#22-database-description)
  * [2.3 Loss Assessment Methods](#23-loss-assessment-methods)
    - [2.3.1 FEMA P-58 Component-Based Method](#231-fema-p-58-component-based-method)
    - [2.3.2 HAZUS Assembly-Based Method](#232-hazus-assembly-based-method)
    - [2.3.3 Story Loss Function (SLF)-Based Method](#233-story-loss-function-slf-based-method)
  * [2.4 Sensitivity Analysis Framework](#24-sensitivity-analysis-framework)
- [3. System Requirements](#3-system-requirements)
- [4. How to Use](#4-how-to-use)
  * [4.1 Comparison of Seismic Loss Assessment Methodologies](#41-comparison-of-seismic-loss-assessment-methodologies)
    - [4.1.1 Steps for Standard Analysis](#411-steps-for-standard-analysis)
  * [4.2 Sensitivity Analysis](#42-sensitivity-analysis)
    - [4.2.1 Piecewise Linear Regression](#421-piecewise-linear-regression)
    - [4.2.2 Weighted EDP Statistics](#422-weighted-edp-statistics)
    - [4.2.3 Residual Drift Sensitivity](#423-residual-drift-sensitivity)
- [5. Repository Contents](#5-repository-contents)
  * [5.1 Scripts](#51-scripts)
  * [5.2 Data](#52-data)
- [6. Credits and Data Sources](#6-credits-and-data-sources)
- [7. Contact](#7-contact)

## 1. Overview and Motivation

Seismic loss estimation is essential tool for assessing building performance and resilience under earthquake hazards. While numerous methodologies exist with varying fidelities, the extent of differences in their estimates and the impact of various modeling decisions remain underexplored. This study addresses these gaps through a systematic comparison of three prominent seismic loss estimation methodologies using a comprehensive database of 621 steel moment-resisting frame buildings.

The research evaluates:
1. **FEMA P-58 component-based approach**: High fidelity and detailed assessment using individual component inventory.
2. **HAZUS assembly-based approach**: Simplified assessment using aggregated assembly inventory  .
3. **Story Loss Function (SLF)-based approach**: The assessment aggregating components at story level. Simplifies using predefined SLFs. Better for portfolio assessments than FEMA P-58

Additionally, sensitivity analyses investigate the influence of key modeling decisions including IM-EDP formulation, EDP proxy selection, demolition fragility characterization, and nonstructural component quantity uncertainties.

## 2. Methodology Overview

### 2.1 Methodology Flowchart

![Methodology Flowchart](https://github.com/shivalingabaddipalli/Loss_Comparison_Sensitivity/blob/main/Figures/Methodology_Flowchart.png)

**Figure 1 (Figure 1 in manuscript). Comparison and sensitivity analysis of seismic loss assessment methodologies.**

The methodology illustrated above shows the comprehensive approach used in this study, encompassing both the comparison of different loss assessment methods and the sensitivity analysis of key modeling decisions. The flowchart demonstrates how building inventory data from the steel building database flows through response analysis, fragility analysis, and loss analysis to ultimately produce Expected Annual Loss (EAL) estimates using different seismic loss estimation methods. Further it presents different sensitivity parameters considered in the study.

### 2.2 Database Description
The analysis utilizes a database of 621 steel special moment resisting frame buildings from Guan et al. (2021) with diverse designs and geometries:

- **Building heights**: 1, 5, 9, 14, and 19 stories
- **Bay configurations**: 1, 3, and 5 moment-resisting frame bays
- **Bay widths**: 20, 30, and 40 feet
- **Site location**: Los Angeles, California (high seismic zone: Ss = 2.25g, S1 = 0.6g)
- **Ground motion suite**: 240 unscaled ground motion records from 12 California earthquakes
- **Structural modeling**: Nonlinear OpenSees models with concentrated plasticity approach

The underlying database is available from Guan et al. (2021) through DesignSafe-CI at https://www.designsafe-ci.org/data/browser/public/designsafe.storage.published/PRJ-2048. 

### 2.3 Loss Assessment Methods

#### 2.3.1 FEMA P-58 Component-Based Method
The most detailed approach evaluating individual building components:
- Structural components: column base plates, splices, beam-column connections
- Non-structural drift-sensitive: partitions, facades, stairs, elevators
- Non-structural acceleration-sensitive: equipment, contents, building systems
- Component quantities determined using FEMA P-58 Normative Tool
- Floor-by-floor analysis using engineering demand parameters

#### 2.3.2 HAZUS Assembly-Based Method  
Simplified approach using aggregated component categories:
- Three main assemblies: structural, non-structural drift-sensitive, non-structural acceleration-sensitive
- Building-level maximum response parameters
- Standardized fragility functions and loss ratios
- Rapid assessment suitable for portfolio evaluation

#### 2.3.3 Story Loss Function (SLF)-Based Method
An approach aggregating component losses at story level:
- Pre-computed loss functions relating EDPs to monetary losses
- Story-level aggregation of component damages
- Reduced computational burden while maintaining component-level detail

### 2.4 Sensitivity Analysis Framework

The study investigates four critical modeling decisions:

1. **IM-EDP Relationship Formulation**: Linear vs. piecewise linear (bilinear) models
2. **EDP Statistics in Assembly-Based Method**: Maximum, mean, weighted averages of percentiles
3. **Demolition Fragility Parameters**: Residual drift threshold and uncertainty characterization
4. **Nonstructural Component Quantities**: Impact of quantity uncertainties on loss estimates

## 3. System Requirements

- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox
- Optimization Toolbox (for weighted EDP statistics and piecewise regression)
- Curve Fitting Toolbox (recommended)

## 4. How to Use

### 4.1 Comparison of Seismic Loss Assessment Methodologies

#### 4.1.1 Steps for Standard Analysis
1. **Prepare workspace**: Ensure all data files are in the working directory and required folders are added to MATLAB path.

2. **Run core loss assessment methods**:
   ```matlab
   % For FEMA P-58 component-based assessment
   FEMAP58_Loss_Assessment.m
   
   % For HAZUS assembly-based assessment  
   HAZUS_Loss_Assessment.m
   
   % For Story Loss Function-based assessment
   SLF_Loss_Assessment.m
   ```

3. **Generate comparison plots**:
   ```matlab
   cd Plotters
   Figure_4_NormalizedEAL.m           % Methodology comparison by building height
   Figure_8_NormalizedEAL_Comparison_EDP_Statistic.m  % EDP statistics comparison
   ```

### 4.2 Sensitivity Analysis

#### 4.2.1 Piecewise Linear Regression

Enhanced IM-EDP demand modeling capturing structural nonlinearity through bilinear relationships.

**Implementation Steps**
1. **Add piecewise functions**: Copy contents of `Piecewise_fit_function.m` to end of main scripts
2. **Replace linear regression calls**: Follow detailed instructions in `Implementation_Guide_Bilinear.md`
3. **Key replacements** (5 locations in each script):
   ```matlab
   % OLD: Linear regression
   md_PIDR = fitlm(log(IM_GM(bd,:)), log(PSDR(flr,:)));
   
   % NEW: Piecewise regression  
   logx_PSDR = log(IM_GM(bd,:))';
   logy_PSDR = log(PSDR(flr,:))';
   best_breakpoint = fminbnd(@(breakpoint) objective(breakpoint, logx_PSDR, logy_PSDR), 
                            min(logx_PSDR), max(logx_PSDR));
   [left_params, right_params, ~, ~, SIGMA_left, SIGMA_right] = 
       piecewise_linear_fit(logx_PSDR, logy_PSDR, best_breakpoint);
   ```

#### 4.2.2 Weighted EDP Statistics

This novel approach improves HAZUS assembly-based estimates by using weighted combinations of EDP percentiles.

**Two-Phase Process**

**Phase 1: Alpha Optimization**
```matlab
% Run optimization to find optimal weighting factors
Alpha_Optimization_EDP_Statistics_25_75_Percentiles.m
```
This script:
- Optimizes weighting factor α for each building
- Uses iterative algorithm with adaptive step size reduction
- Minimizes difference between assembly-based and component-based results
- Outputs: `Alpha_Values_25_75Percentiles.csv`

**Phase 2: Regression Application**  
```matlab
% Apply period-dependent alpha regression model
EDP_Statistics_Sensitivity_analysis.m
```
This script:
- Uses regression equation: `log(α) = f(log(Period))`
- Applies weighted EDP statistics: `Weighted_EDP = α × P75 + (1-α) × P25`
- Computes complete EAL breakdown
- Outputs: Excel file with normalized EAL results

**Supported Modifications for Weighted EDP Analysis**
- **Different percentile combinations**: Modify weighted combination formulas
- **Custom regression models**: Update regression coefficients based on new calibration data
- **Alternative optimization criteria**: Adjust convergence tolerance and step size parameters

#### 4.2.3 Residual Drift Sensitivity

```matlab
% Analyze sensitivity to demolition threshold parameters
RIDR_Sensitivity_analysis.m
```
Evaluates combinations of:
- Demolition thresholds: 0.5% to 5.0% residual drift
- Uncertainty levels: 0.1 to 0.8 logarithmic standard deviation
- Outputs: Grid of results for each parameter combination

## 5. Repository Contents

### 5.1 Scripts
**Core Analysis Scripts:**
- `FEMAP58_Loss_Assessment.m` - FEMA P-58 component-based methodology 
- `HAZUS_Loss_Assessment.m` - HAZUS assembly-based methodology 
- `SLF_Loss_Assessment.m` - Story Loss Function-based methodology 

**Advanced Analysis Scripts:**
- `Alpha_Optimization_EDP_Statistics_25_75_Percentiles.m` - EDP weighting optimization 
- `EDP_Statistics_Sensitivity_analysis.m` - Weighted EDP application and normalized EAL compuatation
- `Piecewise_fit_function.m` - Piecewise regression functions 
- `RIDR_Sensitivity_analysis.m` - Demolition threshold sensitivity 

**Visualization Scripts:**
- `Plotters/Figure_4_NormalizedEAL.m` - Methodology comparison plots
- `Plotters/Figure_8_NormalizedEAL_Comparison_EDP_Statistic.m` - EDP statistics comparison

**Documentation:**
- `Implementation_Guide_Bilinear.md` - Piecewise regression implementation guide 
- `Readme_EDP_Statistic.md` - Weighted EDP statistics documentation 
- `ReadMe_Plotters.txt` - Plotting instructions

### 5.2 Data
**Input Datasets:**
- `GuanDataBase-IMs.xlsx` - Intensity measures for 621 buildings × 240 ground motions 
- `USGSHazard_data.csv` - USGS seismic hazard curves for a site in Los Angeles 
- `Building_Info.xlsx` - Building geometry and period data

**Fragility and Cost Data:**
- `Steel_Building_Data/FEMAP-58_FragilityCostFunctions/` - Component-level fragility and cost functions
- `Steel_Building_Data/HAZUS_FragilityCostFunctions/` - Assembly-level fragility and cost functions  
- `Steel_Building_Data/StoryLossFunctions_20000/` - Pre-computed story loss functions with python codes utilizing the SLF generator library

**Structural Response Data:**
- `Guan_database/StructuralResponses/EDPsUnder240GMs/` - Peak story drift, floor acceleration, and residual drift responses
- `Guan_database/BuildingDesigns/` - Building geometry files

## 6. Credits and Data Sources

- **Steel Building Database**: Guan, M., Burton, H., & Shokrabadi, M. (2021). A database of seismic designs, nonlinear models, and seismic responses for steel moment-resisting frame buildings. *Earthquake Spectra*, 37(2), 1199-1222.
- **FEMA P-58 Methodology**: Federal Emergency Management Agency (2018). Seismic Performance Assessment of Buildings, Volume 1 - Methodology.
- **HAZUS Methodology**: Federal Emergency Management Agency (2015). Multi-hazard Loss Estimation Methodology, Earthquake Model Technical Manual.
- **Seismic Hazard Data**: USGS Unified Hazard Tool (earthquake.usgs.gov/hazards/interactive/)
- **Story Loss Functions**: Based on methodology from Shahnazaryan, D., O'Reilly, G.J., & Monteiro, R. (2021). Story loss functions for seismic design and assessment. *Earthquake Spectra*, 37(4), 2813-2839.

## 7. Contact

**Primary Contact:**
- Shiva Baddipalli (PhD Candidate): shivalinga.baddipalli@usu.edu

**Corresponding Author:**  
- Dr. Mohsen Zaker Esteghamati (Assistant Professor): mohsen.zaker@usu.edu
- Department of Civil and Environmental Engineering
- Utah State University, Logan, UT, United States
