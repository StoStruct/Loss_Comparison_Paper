# Piecewise Linear Regression Implementation Guide

## 1. Objective
Upgrade your main seismic loss assessment linear regression code to use piecewise (bilinear) linear regression for improved accuracy in modeling the relationship between ground motion intensity and structural response. This approach better captures the nonlinear behavior of structures by fitting two linear segments with an optimal breakpoint.

**📋 Note**: While this implementation guide is explained using FEMA P-58 code as the primary example, the same bilinear linear regression approach can be applied to **HAZUS assembly-based** and **Story Loss Function (SLF)-based** methodologies.

## 2. Prerequisites
- The FEMA P-58 linear regression code
- The `piecewise_functions.m` file from this repository
- MATLAB with Optimization Toolbox
  
**📋 Applicable to Multiple Methodologies**: 
- **FEMA P-58**: [`FEMAP58_Loss_Assessment.m`](../FEMAP58_Loss_Assessment.m) (detailed example provided below)
- **HAZUS**: [`HAZUS_Loss_Assessment.m`](../HAZUS_Loss_Assessment.m) (same principles apply)
- **SLF**: [`SLF_Loss_Assessment.m`](../SLF_Loss_Assessment.m) (same principles apply)

## 4. Step-by-Step Implementation

### Step 1: Backup Your Original Code
```matlab
% Create a backup copy of your original script
% For FEMA P-58:
copyfile('../FEMAP58_Loss_Assessment.m', '../FEMAP58_Loss_Assessment_backup.m');

% For HAZUS:
copyfile('../HAZUS_Loss_Assessment.m', '../HAZUS_Loss_Assessment_backup.m');

% For SLF:
copyfile('../SLF_Loss_Assessment.m', '../SLF_Loss_Assessment_backup.m');
```
### Step 2: Add Piecewise Functions
1. Open `Piecewise_fit_function.m` and copy all contents
2. Open your main loss assessment script:
   - FEMA P-58: [`FEMAP58_Loss_Assessment.m`](../FEMAP58_Loss_Assessment.m)
   - HAZUS: [`HAZUS_Loss_Assessment.m`](../HAZUS_Loss_Assessment.m) 
   - SLF: [`SLF_Loss_Assessment.m`](../SLF_Loss_Assessment.m)
3. Paste the piecewise functions **at the very end** of your script, after all main code but before any final `end` statements

⚠️ **Important**: You need to find and replace linear regression calls in your code where `fitlm()` is currently used. The locations and variable names may vary between methodologies:

- **FEMA P-58**: Typically 5 locations (PSDR, PFA, RSDR modeling)
- **HAZUS**: Fewer locations due to simplified approach (typically 2-3 locations)
- **SLF**: Variable depending on implementation (typically 3-4 locations)

**📋 Implementation Strategy**: The examples below show FEMA P-58 implementation. For HAZUS and SLF codes:
1. Look for similar `fitlm()` patterns in your code
2. Apply the same replacement logic
3. Adjust variable names to match your specific code structure

For FEMA P-58 code, you need to find and replace **5 different locations** in your code where linear regression is used:

#### Location 1: Structural Components (PSDR modeling)
**Find this pattern around line 190-200:**
```matlab
md_PIDR = fitlm(log(IM_GM(bd, :)), log(PSDR(flr, :)));
PIDR_DemandPar(1:2) = md_PIDR.Coefficients{1:2, 1};
PIDR_DemandPar(3) = md_PIDR.RMSE;
SIGMA = PIDR_DemandPar(3);
MD_PIDR = PIDR_DemandPar(2) * log(SaCalc(im)) + PIDR_DemandPar(1);
```
**Replace with:**
```matlab
% Piecewise linear regression for PSDR
logx_PSDR = log(IM_GM(bd, :))';
logy_PSDR = log(PSDR(flr,:))';
best_breakpoint_PSDR = fminbnd(@(breakpoint) objective(breakpoint, logx_PSDR, logy_PSDR), ...
                               min(logx_PSDR), max(logx_PSDR));
[left_params_PSDR, right_params_PSDR, ~, ~, SIGMA_left_PSDR, SIGMA_right_PSDR] = ...
    piecewise_linear_fit(logx_PSDR, logy_PSDR, best_breakpoint_PSDR);

% Select appropriate segment based on intensity level
if log(SaCalc(im)) <= best_breakpoint_PSDR
    MD_PIDR = left_params_PSDR(1) * log(SaCalc(im)) + left_params_PSDR(2);
    SIGMA = SIGMA_left_PSDR;
else
    MD_PIDR = right_params_PSDR(1) * log(SaCalc(im)) + right_params_PSDR(2);
    SIGMA = SIGMA_right_PSDR;
end
```
#### Location 2: Non-Structural Drift-Sensitive Components 
**Find similar pattern around line 300-310 in the drift-sensitive components section**

Apply the same replacement logic as Location 1, ensuring variable names match the context.

#### Location 3: Non-Structural Acceleration-Sensitive (Floor Level)
**Find this pattern around line 400-410:**

```matlab
md_PFA = fitlm(log(IM_GM(bd, :)), log(PFA(flr, :)));
PFA_DemandPar(1:2) = md_PFA.Coefficients{1:2, 1};
PFA_DemandPar(3) = md_PFA.RMSE;
SIGMA_PFA = PFA_DemandPar(3);
MD_PFA = PFA_DemandPar(2) * log(SaCalc(im)) + PFA_DemandPar(1);
```

**Replace with:**
```matlab
% Piecewise linear regression for PFA
logx_PFA = log(IM_GM(bd, :))';
logy_PFA = log(PFA(flr, :))';
best_breakpoint_PFA = fminbnd(@(breakpoint) objective(breakpoint, logx_PFA, logy_PFA), ...
                              min(logx_PFA), max(logx_PFA));
[left_params_PFA, right_params_PFA, ~, ~, SIGMA_left_PFA, SIGMA_right_PFA] = ...
    piecewise_linear_fit(logx_PFA, logy_PFA, best_breakpoint_PFA);

% Select appropriate segment
if log(SaCalc(im)) <= best_breakpoint_PFA
    MD_PFA = left_params_PFA(1) * log(SaCalc(im)) + left_params_PFA(2);
    SIGMA_PFA = SIGMA_left_PFA;
else
    MD_PFA = right_params_PFA(1) * log(SaCalc(im)) + right_params_PFA(2);
    SIGMA_PFA = SIGMA_right_PFA;
end
```


#### Location 4: Building-Level Acceleration Components
**Find similar PFA modeling pattern around line 500-510**

Apply the same PFA replacement as Location 3.

#### Location 5: Demolition Loss (RSDR modeling)
**Find this pattern around line 600-610:**
```matlab
md_RSDR = fitlm(log(IM_GM(bd,:)), log(Max_RSDR(1,:)));
RSDR_DemandPar(1:2) = md_RSDR.Coefficients{1:2, 1};
RSDR_DemandPar(3) = md_RSDR.RMSE;
SIGMA_RSDR = RSDR_DemandPar(3);
MD_RSDR = RSDR_DemandPar(2)*log(SaCalc(im))+RSDR_DemandPar(1);
```

**Replace with:**
```matlab
% Piecewise linear regression for RSDR
logx_RSDR = log(IM_GM(bd,:))';
logy_RSDR = log(Max_RSDR(1,:))';
best_breakpoint_RSDR = fminbnd(@(breakpoint) objective(breakpoint, logx_RSDR, logy_RSDR), ...
                               min(logx_RSDR), max(logx_RSDR));
[left_params_RSDR, right_params_RSDR, ~, ~, SIGMA_left_RSDR, SIGMA_right_RSDR] = ...
    piecewise_linear_fit(logx_RSDR, logy_RSDR, best_breakpoint_RSDR);

% Select appropriate segment
if log(SaCalc(im)) <= best_breakpoint_RSDR
    MD_RSDR = left_params_RSDR(1) * log(SaCalc(im)) + left_params_RSDR(2);
    SIGMA_RSDR = SIGMA_left_RSDR;
else
    MD_RSDR = right_params_RSDR(1) * log(SaCalc(im)) + right_params_RSDR(2);
    SIGMA_RSDR = SIGMA_right_RSDR;
end
```

