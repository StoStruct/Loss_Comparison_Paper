%% PIECEWISE LINEAR REGRESSION FUNCTIONS FOR FEMA P-58
% =========================================================================
% These functions enable piecewise (bilinear) linear regression for 
% improved demand modeling in FEMA P-58 seismic loss assessment
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
% HOW TO USE:
% 1. Add these functions to the end of your main FEMA P-58 script
% 2. Follow the modification guide in IMPLEMENTATION_GUIDE.md
% 3. Replace all fitlm() calls with piecewise regression blocks
% =========================================================================

%% Function 1: Piecewise Linear Fit
function [left_params, right_params, left_eq, right_eq, SIGMA_left, SIGMA_right] = piecewise_linear_fit(logx, logy, breakpoint)
    % PIECEWISE_LINEAR_FIT Fits two linear segments given a breakpoint
    %
    % INPUTS:
    %   logx - Log of independent variable (e.g., log(Sa))
    %   logy - Log of dependent variable (e.g., log(PSDR), log(PFA))
    %   breakpoint - Optimal breakpoint location
    %
    % OUTPUTS:
    %   left_params - [slope, intercept] for left segment
    %   right_params - [slope, intercept] for right segment  
    %   left_eq - String equation for left segment
    %   right_eq - String equation for right segment
    %   SIGMA_left - Standard deviation for left segment
    %   SIGMA_right - Standard deviation for right segment

    % Split the data into two segments based on the breakpoint
    left_mask = logx <= breakpoint;
    right_mask = logx > breakpoint;

    % Linear fit for the left segment
    left_x = logx(left_mask);
    left_y = logy(left_mask);
    left_params = polyfit(left_x, left_y, 1);

    % Linear fit for the right segment
    right_x = logx(right_mask);
    right_y = logy(right_mask);
    right_params = polyfit(right_x, right_y, 1);

    % Ensure continuity at the breakpoint
    y_breakpoint = polyval(left_params, breakpoint);
    right_params(2) = y_breakpoint - right_params(1) * breakpoint;

    % Calculate SIGMA for both segments
    residuals_left = left_y - polyval(left_params, left_x);
    residuals_right = right_y - polyval(right_params, right_x);
    SIGMA_left = std(residuals_left);
    SIGMA_right = std(residuals_right);

    % Generate equations for the left and right segments (optional)
    left_eq = sprintf('log(y) = %.4f * log(x) + %.4f', left_params(1), left_params(2));
    right_eq = sprintf('log(y) = %.4f * log(x) + %.4f', right_params(1), right_params(2));
end

%% Function 2: Objective Function for Optimization
function ssr = objective(breakpoint, logx, logy)
    % OBJECTIVE Objective function to minimize for finding optimal breakpoint
    %
    % INPUTS:
    %   breakpoint - Candidate breakpoint location
    %   logx - Log of independent variable
    %   logy - Log of dependent variable
    %
    % OUTPUT:
    %   ssr - Sum of squared residuals for the piecewise fit

    % Get the fitted parameters for this breakpoint
    [left_params, right_params, ~, ~, ~, ~] = piecewise_linear_fit(logx, logy, breakpoint);

    % Calculate residuals for both segments
    left_mask = logx <= breakpoint;
    right_mask = logx > breakpoint;

    % Compute residuals for both segments
    left_residuals = logy(left_mask) - polyval(left_params, logx(left_mask));
    right_residuals = logy(right_mask) - polyval(right_params, logx(right_mask));

    % Combine all residuals
    total_residuals = [left_residuals; right_residuals];

    % Calculate sum of squared residuals (objective to minimize)
    ssr = sum(total_residuals .^ 2);
end

%% USAGE EXAMPLE:
% Replace this old block with new block in the main code 
% (FEMAP58_Loss_Assesment) to implement bilinear or Piecewise Linear Regression:
%
% OLD CODE (Linear Regression in FEMAP58_Loss_Assesment) :
% md_PIDR = fitlm(log(IM_GM(bd, :)), log(PSDR(flr, :)));
% PIDR_DemandPar(1:2) = md_PIDR.Coefficients{1:2, 1};
% PIDR_DemandPar(3) = md_PIDR.RMSE;
% SIGMA = PIDR_DemandPar(3);
% MD_PIDR = PIDR_DemandPar(2) * log(SaCalc(im)) + PIDR_DemandPar(1);
%
% NEW CODE (Piecewise Linear Regression):
% logx_PSDR = log(IM_GM(bd, :))';
% logy_PSDR = log(PSDR(flr,:))';
% best_breakpoint_PSDR = fminbnd(@(breakpoint) objective(breakpoint, logx_PSDR, logy_PSDR), min(logx_PSDR), max(logx_PSDR));
% [left_params_PSDR, right_params_PSDR, ~, ~, SIGMA_left_PSDR, SIGMA_right_PSDR] = piecewise_linear_fit(logx_PSDR, logy_PSDR, best_breakpoint_PSDR);
% 
% if log(SaCalc(im)) <= best_breakpoint_PSDR
%     MD_PIDR = left_params_PSDR(1) * log(SaCalc(im)) + left_params_PSDR(2);
%     SIGMA = SIGMA_left_PSDR;
% else
%     MD_PIDR = right_params_PSDR(1) * log(SaCalc(im)) + right_params_PSDR(2);
%     SIGMA = SIGMA_right_PSDR;
% end
