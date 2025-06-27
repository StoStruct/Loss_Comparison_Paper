# Piecewise Linear Regression Implementation Guide

## Objective
Upgrade your FEMA P-58 linear regression code to use piecewise (bilinear) linear regression for improved accuracy in modeling the relationship between ground motion intensity and structural response.

## Prerequisites
- The FEMA P-58 linear regression code
- The `piecewise_functions.m` file from this repository
- MATLAB with Optimization Toolbox

## Step-by-Step Implementation

### Step 1: Add Piecewise Functions
Copy the entire contents of `piecewise_functions.m` and paste it **at the very end** of your main FEMA P-58 script, after all the main code but before the final `end` statement.

### Step 2: Locate and Replace Linear Regression Calls

You need to find and replace **5 different locations** in your code where linear regression is used:

#### Location 1: Structural Components (PSDR modeling)
**Find this pattern around line 190-200:**

md_PIDR = fitlm(log(IM_GM(bd, :)), log(PSDR(flr, :)));
PIDR_DemandPar(1:2) = md_PIDR.Coefficients{1:2, 1};
PIDR_DemandPar(3) = md_PIDR.RMSE;
SIGMA = PIDR_DemandPar(3);
MD_PIDR = PIDR_DemandPar(2) * log(SaCalc(im)) + PIDR_DemandPar(1);


**Replace with:**

% Piecewise linear regression for PSDR
logx_PSDR = log(IM_GM(bd, :))';
logy_PSDR = log(PSDR(flr,:))';
best_breakpoint_PSDR = fminbnd(@(breakpoint) objective(breakpoint, logx_PSDR, logy_PSDR), min(logx_PSDR), max(logx_PSDR));
[left_params_PSDR, right_params_PSDR, ~, ~, SIGMA_left_PSDR, SIGMA_right_PSDR] = piecewise_linear_fit(logx_PSDR, logy_PSDR, best_breakpoint_PSDR);

% Select appropriate segment based on intensity level
if log(SaCalc(im)) <= best_breakpoint_PSDR
    MD_PIDR = left_params_PSDR(1) * log(SaCalc(im)) + left_params_PSDR(2);
    SIGMA = SIGMA_left_PSDR;
else
    MD_PIDR = right_params_PSDR(1) * log(SaCalc(im)) + right_params_PSDR(2);
    SIGMA = SIGMA_right_PSDR;
end

#### Location 2: Non-Structural Drift-Sensitive Components 
**Find similar pattern around line 300-310 in the drift-sensitive components section**

Apply the same replacement as Location 1.

#### Location 3: Non-Structural Acceleration-Sensitive (Floor Level)
**Find this pattern around line 400-410:**

md_PFA = fitlm(log(IM_GM(bd, :)), log(PFA(flr, :)));
PFA_DemandPar(1:2) = md_PFA.Coefficients{1:2, 1};
PFA_DemandPar(3) = md_PFA.RMSE;
SIGMA_PFA = PFA_DemandPar(3);
MD_PFA = PFA_DemandPar(2) * log(SaCalc(im)) + PFA_DemandPar(1);

**Replace with:**

% Piecewise linear regression for PFA
logx_PFA = log(IM_GM(bd, :))';
logy_PFA = log(PFA(flr, :))';
best_breakpoint_PFA = fminbnd(@(breakpoint) objective(breakpoint, logx_PFA, logy_PFA), min(logx_PFA), max(logx_PFA));
[left_params_PFA, right_params_PFA, ~, ~, SIGMA_left_PFA, SIGMA_right_PFA] = piecewise_linear_fit(logx_PFA, logy_PFA, best_breakpoint_PFA);

% Select appropriate segment
if log(SaCalc(im)) <= best_breakpoint_PFA
    MD_PFA = left_params_PFA(1) * log(SaCalc(im)) + left_params_PFA(2);
    SIGMA_PFA = SIGMA_left_PFA;
else
    MD_PFA = right_params_PFA(1) * log(SaCalc(im)) + right_params_PFA(2);
    SIGMA_PFA = SIGMA_right_PFA;
end


#### Location 4: Building-Level Acceleration Components
**Find similar PFA modeling pattern around line 500-510**

Apply the same PFA replacement as Location 3.

#### Location 5: Demolition Loss (RSDR modeling)
**Find this pattern around line 600-610:**

md_RSDR = fitlm(log(IM_GM(bd,:)), log(Max_RSDR(1,:)));
RSDR_DemandPar(1:2) = md_RSDR.Coefficients{1:2, 1};
RSDR_DemandPar(3) = md_RSDR.RMSE;
SIGMA_RSDR = RSDR_DemandPar(3);
MD_RSDR = RSDR_DemandPar(2)*log(SaCalc(im))+RSDR_DemandPar(1);


**Replace with:**

% Piecewise linear regression for RSDR
logx_RSDR = log(IM_GM(bd,:))';
logy_RSDR = log(Max_RSDR(1,:))';
best_breakpoint_RSDR = fminbnd(@(breakpoint) objective(breakpoint, logx_RSDR, logy_RSDR), min(logx_RSDR), max(logx_RSDR));
[left_params_RSDR, right_params_RSDR, ~, ~, SIGMA_left_RSDR, SIGMA_right_RSDR] = piecewise_linear_fit(logx_RSDR, logy_RSDR, best_breakpoint_RSDR);

% Select appropriate segment
if log(SaCalc(im)) <= best_breakpoint_RSDR
    MD_RSDR = left_params_RSDR(1) * log(SaCalc(im)) + left_params_RSDR(2);
    SIGMA_RSDR = SIGMA_left_RSDR;
else
    MD_RSDR = right_params_RSDR(1) * log(SaCalc(im)) + right_params_RSDR(2);
    SIGMA_RSDR = SIGMA_right_RSDR;
end

