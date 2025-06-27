clc;
clear;
close all;

% Load the data
file1 = 'MyInventory_Structural_4comp\CombinedOutput_SLF_Struct_4comps.csv';
% file1 = 'MyInventory_Structural_4comp - Copy\CombinedOutput_SLF_Struct_4comps.csv'; %% DONT USE THIS FOLDER BECAUSE SLAB REPAIR COST FUNCTIONS WERE WRONG
fourcomp = readmatrix(file1);

% Extract columns for IDR and SLFs
IDR = fourcomp(:,1);
SLF_1 = fourcomp(:,2);
SLF_2 = fourcomp(:,3);
SLF_3 = fourcomp(:,4);
SLF_4 = fourcomp(:,5);
SLF_5 = fourcomp(:,6);
SLF_6 = fourcomp(:,7);
SLF_7 = fourcomp(:,8);
SLF_8 = fourcomp(:,9);
SLF_9 = fourcomp(:,10);

% Set the font to Times New Roman
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');

% Create figure
figure('Position', [100, 100, 400, 350]);
hold on;

% % Plotting each SLF against IDR
% plot(IDR, SLF_1, 'LineWidth', 2, 'DisplayName', 'SLF 1');
% plot(IDR, SLF_2, 'LineWidth', 2, 'DisplayName', 'SLF 2', 'LineStyle', '--');
% plot(IDR, SLF_3, 'LineWidth', 2, 'DisplayName', 'SLF 3', 'LineStyle', ':');
% plot(IDR, SLF_4, 'LineWidth', 2, 'DisplayName', 'SLF 4', 'LineStyle', '-.');
% plot(IDR, SLF_5, 'LineWidth', 2, 'DisplayName', 'SLF 5', 'LineStyle', '-');
% plot(IDR, SLF_6, 'LineWidth', 2, 'DisplayName', 'SLF 6', 'LineStyle', '--');
% plot(IDR, SLF_7, 'LineWidth', 2, 'DisplayName', 'SLF 7', 'LineStyle', ':');
% plot(IDR, SLF_8, 'LineWidth', 2, 'DisplayName', 'SLF 8', 'LineStyle', '-.');
% plot(IDR, SLF_9, 'LineWidth', 2, 'DisplayName', 'SLF 9', 'LineStyle', '-');


% % Plotting each SLF against IDR
% plot(IDR, SLF_1, 'LineWidth', 2, 'DisplayName', 'BC-1, BW-20');
% % plot(IDR, SLF_2, 'LineWidth', 2, 'DisplayName', 'BC-1, BW-30', 'LineStyle', '--');
% plot(IDR, SLF_3, 'LineWidth', 2, 'DisplayName', 'BC-1, BW-40', 'LineStyle', ':');
% plot(IDR, SLF_4, 'LineWidth', 2, 'DisplayName', 'BC-3, BW-20', 'LineStyle', '-.');
% % plot(IDR, SLF_5, 'LineWidth', 2, 'DisplayName', 'BC-3, BW-30', 'LineStyle', '-');
% plot(IDR, SLF_6, 'LineWidth', 2, 'DisplayName', 'BC-3, BW-40', 'LineStyle', '--');
% plot(IDR, SLF_7, 'LineWidth', 2, 'DisplayName', 'BC-5, BW-20', 'LineStyle', ':');
% % plot(IDR, SLF_8, 'LineWidth', 2, 'DisplayName', 'BC-5, BW-30', 'LineStyle', '-.');
% plot(IDR, SLF_9, 'LineWidth', 2, 'DisplayName', 'BC-5, BW-40', 'LineStyle', '-');

% Define a custom colormap with enough distinct colors
colorSet = [
    0.1216, 0.4667, 0.7059;  % blue
    1.0000, 0.4980, 0.0549;  % orange
    0.1725, 0.6275, 0.1725;  % green
    0.8392, 0.1529, 0.1569;  % red
    0.5804, 0.4039, 0.7412;  % purple
    0.5490, 0.3373, 0.2941;  % brown
    0.8902, 0.4667, 0.7608   % pink
];

hold on

% Plot each SLF with a unique color
plot(IDR, SLF_1, 'LineWidth', 2, 'Color', colorSet(1,:), 'DisplayName', 'BC-1, BW-20');
plot(IDR, SLF_3, 'LineWidth', 2, 'LineStyle', ':',  'Color', colorSet(2,:), 'DisplayName', 'BC-1, BW-40');
plot(IDR, SLF_4, 'LineWidth', 2, 'LineStyle', '-.', 'Color', colorSet(3,:), 'DisplayName', 'BC-3, BW-20');
plot(IDR, SLF_6, 'LineWidth', 2, 'LineStyle', '--', 'Color', colorSet(4,:), 'DisplayName', 'BC-3, BW-40');
plot(IDR, SLF_7, 'LineWidth', 2, 'LineStyle', ':',  'Color', colorSet(5,:), 'DisplayName', 'BC-5, BW-20');
plot(IDR, SLF_9, 'LineWidth', 2, 'LineStyle', '-',  'Color', colorSet(6,:), 'DisplayName', 'BC-5, BW-40');

% legend('Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 10);

% Adding labels and title
xlabel('PIDR', 'FontSize', 19, 'Color', 'k');
ylabel('Story loss (USD)', 'FontSize', 19, 'Color', 'k');
% title('SLFs for Column Base Plates', 'FontSize', 14);

% % Enhance the plot box
% set(gca, 'LineWidth', 1.5); % Make the box thicker
% set(gca, 'FontSize', 12); % Setting font size for axes
% set(gca, 'Box', 'on'); % Ensure the box is on
% set(gca, 'LineWidth', 1.2, 'Box', 'on', 'XColor', 'k', 'YColor', 'k');
% Adjust tick labels font size and axes properties
set(gca, 'FontSize', 19, 'LineWidth', 1.5, 'Box', 'on', 'XColor', 'k', 'YColor', 'k');

legend('show', 'Location', 'northeastoutside', 'FontSize', 16);

% Grid on for better readability
% grid on;

% Set the axes for better visualization

ylim([0 2000000]); % Adjust if necessary
xlim([0 0.2]);
xticks([0, 0.05, 0.1, 0.15, 0.2]);  % Set y-axis ticks

% Save the figure

print(gcf, 'SLF_vs_IDR_4comp_legend', '-dpng', '-r1200');


% clc;
% clear;
% close all;
% 
% % Load the data
% file1 = 'MyInventory_Structural_4comp\CombinedOutput_SLF_Struct_4comps.csv';
% % file1 = 'MyInventory_Structural_4comp - Copy\CombinedOutput_SLF_Struct_4comps.csv'; %% DONT USE THIS FOLDER BECAUSE SLAB REPAIR COST FUNCTIONS WERE WRONG
% fourcomp = readmatrix(file1);
% 
% % Extract columns for IDR and SLFs
% IDR = fourcomp(:,1);
% SLF_1 = fourcomp(:,2);
% SLF_2 = fourcomp(:,3);
% SLF_3 = fourcomp(:,4);
% SLF_4 = fourcomp(:,5);
% SLF_5 = fourcomp(:,6);
% SLF_6 = fourcomp(:,7);
% SLF_7 = fourcomp(:,8);
% SLF_8 = fourcomp(:,9);
% SLF_9 = fourcomp(:,10);
% 
% % Set the font to Times New Roman
% set(0, 'DefaultAxesFontName', 'Times New Roman');
% set(0, 'DefaultTextFontName', 'Times New Roman');
% 
% % Create figure
% figure('Position', [100, 100, 440, 350]);
% hold on;
% 
% % Plotting each SLF against IDR
% plot(IDR, SLF_1, 'LineWidth', 2, 'DisplayName', 'SLF 1');
% plot(IDR, SLF_2, 'LineWidth', 2, 'DisplayName', 'SLF 2', 'LineStyle', '--');
% plot(IDR, SLF_3, 'LineWidth', 2, 'DisplayName', 'SLF 3', 'LineStyle', ':');
% plot(IDR, SLF_4, 'LineWidth', 2, 'DisplayName', 'SLF 4', 'LineStyle', '-.');
% plot(IDR, SLF_5, 'LineWidth', 2, 'DisplayName', 'SLF 5', 'LineStyle', '-');
% plot(IDR, SLF_6, 'LineWidth', 2, 'DisplayName', 'SLF 6', 'LineStyle', '--');
% plot(IDR, SLF_7, 'LineWidth', 2, 'DisplayName', 'SLF 7', 'LineStyle', ':');
% plot(IDR, SLF_8, 'LineWidth', 2, 'DisplayName', 'SLF 8', 'LineStyle', '-.');
% plot(IDR, SLF_9, 'LineWidth', 2, 'DisplayName', 'SLF 9', 'LineStyle', '-');
% 
% % Adding labels and title
% xlabel('PIDR', 'FontSize', 19, 'Color', 'k');
% ylabel('Story loss function (USD)', 'FontSize', 19, 'Color', 'k');
% % title('SLFs for Column Base Plates', 'FontSize', 14);
% 
% % Enhance the plot box
% set(gca, 'FontSize', 19, 'LineWidth', 1.5, 'Box', 'on', 'XColor', 'k', 'YColor', 'k');
% 
% % Place the legend horizontally in two lines
% legend_handle = legend('show', 'Orientation', 'horizontal', 'FontSize', 14);
% legend_handle.NumColumns = 5; % Stack horizontally in two lines with 5 items per row
% set(legend_handle, 'Location', 'north');
% % Grid on for better readability
% % grid on;
% 
% % Set the axes for better visualization
% axis tight;
% ylim([0 2000000]); % Adjust if necessary
% xlim([0 0.2]);
% 
% % Save the figure
% saveas(gcf, 'SLF_vs_IDR_4comp.png'); % Save as PNG
% 
