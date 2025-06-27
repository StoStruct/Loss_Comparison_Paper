%% ========================================================================
%  FIGURE 8A: COMPARISON OF DIFFERENT EDP statistic FOR TOTAL EAL
%  
%  This script creates box plots comparing different HAZUS assembly-based
%  EDP statistic against component-based and SLF approaches for total EAL.
%
%  Methodologies compared:
%  1. Assembly (Max) - HAZUS maximum EDP approach
%  2. Assembly (Mean) - HAZUS mean EDP approach  
%  3. Assembly (WA_Max-Mean) - Weighted average of max and mean
%  4. Assembly (WA_25-75) - Weighted average of 25th and 75th percentiles
%  5. Component - FEMA P-58 component-based approach
%  6. SLF - Story Loss Function approach
%
%  Data Column: Column 11 (Total losses)
%
% % Author: Shiva Baddipalli
% Institution: Utah State University, USA
% Email: shivalinga.baddipalli@usu.edu
% Last Updated: June 25, 2025
%% ========================================================================

%% Clear workspace and close all figures
clc;
clearvars; 
close all;
clear;

%% File paths - Update these paths according to your data location
csvFilePath1 = 'Building_Info.xlsx';
Hazus_Max_FILE = 'HAZUS_EAL_NormEAL.xlsx';
Hazus_Mean_FILE = 'HAZUS_NormEAL_Mean.xlsx';
Hazus_WA_FILE_TEST_MaxMean='HAZUS_NormEAL_WeightedAverage_MaxMean.xlsx'; 
Hazus_WA_FILE_TEST_25_75='HAZUS_NormEAL_WeightedAverage_25-75Percentile.xlsx';
csvFilePath2 = 'FEMAP58_EAL_NormEAL.xlsx';
SLF_FILE = 'SLF_EAL_NormEAL.xlsx';

%% Read data from Excel files
% Building information data
data1 = readmatrix(csvFilePath1);
data1 = data1(2:end, :);                    % Remove header row

% EAL data from different methodologies
data2 = readmatrix(csvFilePath2);           % FEMA P-58 component-based data
data2_hazus_max = readmatrix(Hazus_Max_FILE);       % HAZUS maximum EDP data
data_hazus_mean = readmatrix(Hazus_Mean_FILE);      % HAZUS mean EDP data

data_hazus_WAT_MaxMean = readmatrix(Hazus_WA_FILE_TEST_MaxMean);    % Weighted average: Max-Mean

data_hazus_WAT_25_75 = readmatrix(Hazus_WA_FILE_TEST_25_75);       % Weighted average: 25th-75th percentile
data_SLF = readmatrix(SLF_FILE);            % Story Loss Function data

%% Extract normalized EAL data from column 11 (Total losses) and convert to percentage
% NOTE: All data uses column 11 for total losses
NormEAL_Tot_hazus_max = data2_hazus_max(:, 11).*100;        % HAZUS Maximum approach
NormEAL_Tot = data2(:, 11).*100;                            % Component-based (FEMA P-58)
NormEAL_Tot_hazus_mean = data_hazus_mean(:, 11).*100;       % HAZUS Mean approach

NormEAL_Tot_hazus_WAT_MaxMean = data_hazus_WAT_MaxMean(:, 11).*100;  % Weighted average: Max-Mean

NormEAL_Tot_hazus_WAT_25_75 = data_hazus_WAT_25_75(:, 11).*100;     % Weighted average: 25th-75th percentile
NormEAL_Tot_SLF = data_SLF(:, 11).*100;                     % Story Loss Function

%% Define color scheme for different methodologies
% Primary methodology colors
dark_red = [0.937, 0.231, 0.173];      % Component-based (FEMA P-58)
dark_blue = [0.259, 0.573, 0.776];     % HAZUS Maximum approach
dark_green = [0.451, 0.451, 0.451];    % Story Loss Function

% Gradient blues for different HAZUS approaches
light_blue = [0.42, 0.682, 0.839];           % HAZUS Mean approach
lightLight_blue = [0.62, 0.792, 0.882];      % HAZUS Weighted Average Max-Mean
lightLightLight_blue = [0.776, 0.859, 0.937]; % HAZUS Weighted Average 25th-75th

%% Create figure with specific dimensions for publication
figure('Position', [100, 100, 510, 355]);
hold on;

%% Create box plots for each methodology (positioned from 1 to 6)

% Position 1: Assembly-based (Maximum EDP approach)
h1 = boxplot(NormEAL_Tot_hazus_max, 'Colors', dark_blue, 'Positions', 1, 'Widths', 0.7, 'Symbol', '+');
patch(get(findobj(gca,'Tag','Box','-and','Color',dark_blue),'XData'), ...
      get(findobj(gca,'Tag','Box','-and','Color',dark_blue),'YData'), dark_blue, 'FaceAlpha', 0.7);
hold on;

% Position 2: Assembly-based (Mean EDP approach)
h2 = boxplot(NormEAL_Tot_hazus_mean, 'Colors', light_blue, 'Positions', 2, 'Widths', 0.7, 'Symbol', '+');
patch(get(findobj(gca,'Tag','Box','-and','Color',light_blue),'XData'), ...
      get(findobj(gca,'Tag','Box','-and','Color',light_blue),'YData'), light_blue, 'FaceAlpha', 0.7);
hold on;

% Position 3: Assembly-based (Weighted Average - Max and Mean)
h3 = boxplot(NormEAL_Tot_hazus_WAT_MaxMean, 'Colors', lightLight_blue, 'Positions', 3, 'Widths', 0.7, 'Symbol', '+');
patch(get(findobj(gca,'Tag','Box','-and','Color',lightLight_blue),'XData'), ...
      get(findobj(gca,'Tag','Box','-and','Color',lightLight_blue),'YData'), lightLight_blue, 'FaceAlpha', 0.7);
hold on;

% Position 4: Assembly-based (Weighted Average - 25th and 75th Percentile)
h4 = boxplot(NormEAL_Tot_hazus_WAT_25_75, 'Colors', lightLightLight_blue, 'Positions', 4, 'Widths', 0.7, 'Symbol', '+');
patch(get(findobj(gca,'Tag','Box','-and','Color',lightLightLight_blue),'XData'), ...
      get(findobj(gca,'Tag','Box','-and','Color',lightLightLight_blue),'YData'), lightLightLight_blue, 'FaceAlpha', 0.7);
hold on;

% Position 5: Component-based (FEMA P-58)
h5 = boxplot(NormEAL_Tot, 'Colors', dark_red, 'Positions', 5, 'Widths', 0.7, 'Symbol', '+');
patch(get(findobj(gca,'Tag','Box','-and','Color',dark_red),'XData'), ...
      get(findobj(gca,'Tag','Box','-and','Color',dark_red),'YData'), dark_red, 'FaceAlpha', 0.7);
hold on;

% Position 6: SLF-based (Story Loss Function)
h6 = boxplot(NormEAL_Tot_SLF, 'Colors', dark_green, 'Positions', 6, 'Widths', 0.7, 'Symbol', '+');
patch(get(findobj(gca,'Tag','Box','-and','Color',dark_green),'XData'), ...
      get(findobj(gca,'Tag','Box','-and','Color',dark_green),'YData'), dark_green, 'FaceAlpha', 0.7);

%% Set axis limits and ticks for total EAL comparison
ylim([0, 15]);          % Y-axis range for total EAL percentage
xlim([0.5, 6.5]);       % X-axis range to accommodate 6 methodologies
yticks([0, 5, 10, 15]); % Y-axis tick marks

%% Add axis labels
ylabel('Normalized EAL (%)', 'FontSize', 18, 'Color', 'k');

%% Set x-axis tick labels with methodology names
% Option 1: Full descriptive labels (currently commented out)
% set(gca, 'XTick', [1, 2, 3, 4, 5, 6], 'XTickLabel', ...
%     {'Assembly (Max)', 'Assembly (Mean)', 'Assembly (WA - Max and Mean)', ...
%     'Assembly (WA - 25^{th} and 75^{th} Percentile)', 'Component', 'SLF'}, ...
%     'TickLabelInterpreter', 'tex', 'FontSize', 20);

% Option 2: Labels with subscripts (currently commented out)
% set(gca, 'XTick', [1, 2, 3, 4, 5, 6], 'XTickLabel', ...
%     {'Assembly (Max)', 'Assembly (Mean)', 'Assembly (WA_{Max-Mean})', ...
%     'Assembly (WA_{25^{th}-75})', 'Component', 'SLF'}, ...
%     'TickLabelInterpreter', 'tex', 'FontSize', 20);

% Option 3: Compact labels with subscripts (currently active)
set(gca, 'XTick', [1, 2, 3, 4, 5, 6], 'XTickLabel', ...
    {'Max', 'Mean', 'WA_{Max-Mean}', ...
    'WA_{25-75}', 'Component', 'SLF'}, ...
    'TickLabelInterpreter', 'tex', 'FontSize', 18);

%% Adjust line thickness for all box plot elements
set(findobj(gca,'type','line'),'linew',1.75)

%% Enhance median lines appearance
% Bring median lines to front and make them black with increased thickness
medians = findobj(gca, 'Tag', 'Median');
for i = 1:length(medians)
    uistack(medians(i), 'top');                    % Bring median line to top
    set(medians(i), 'Color', 'k', 'LineWidth', 1.75); % Set color to black and increase thickness
end

%% Format axes with professional appearance
set(gca, 'LineWidth', 1.5, 'Box', 'on', 'XColor', 'k', 'YColor', 'k');

%% Add horizontal grid lines for better readability
grid on;
set(gca, 'XGrid', 'off', 'YGrid', 'on', 'GridLineStyle', '--', ...
         'GridColor', [0.8, 0.8, 0.8], 'GridAlpha', 1);

%% Optional: Additional grid formatting (currently commented out)
% grid on;
% set(gca, 'XGrid', 'off', 'YGrid', 'on', 'GridLineStyle', '--', 'GridColor', [0.8, 0.8, 0.8]);

%% Save figure options (currently commented out)
% Update filename as needed
% saveas(gcf, 'MaxvsMeanvsWAT_MaxMeanvsWAT_25_75_Filled_10-19-24.png');
% print(gcf, 'Figure_8a_04-12-25.jpg', '-djpeg', '-r2400');