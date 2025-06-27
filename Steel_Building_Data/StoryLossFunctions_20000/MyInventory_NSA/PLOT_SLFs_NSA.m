clc;
clear;
close all;

% Load the data
file1 = 'MyInventory_NSA_FloorLevelComponents\CombinedOutput_SLF_NSA_Flrlvl.csv';
NSA = readmatrix(file1);

% Extract columns for IDR and SLFs
IDR = NSA(:,1);
SLF_1 = NSA(:,2);
SLF_2 = NSA(:,3);
SLF_3 = NSA(:,4);

% Set the font to Times New Roman
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');

% Create figure
figure('Position', [100, 100, 390, 350]);
hold on;

% Plotting each SLF against IDR
plot(IDR, SLF_1, 'LineWidth', 2, 'DisplayName', 'BW-20');
plot(IDR, SLF_2, 'LineWidth', 2, 'DisplayName', 'BW-30', 'LineStyle', '--');
plot(IDR, SLF_3, 'LineWidth', 2, 'DisplayName', 'BW-40', 'LineStyle', ':');

% Adding labels and title
% xlabel('PFA (g)', 'FontSize', 19, 'Color', 'k');
% xlabel('PFA (\it{g})', 'FontSize', 19, 'Color', 'k', 'Interpreter', 'latex');
xlabel('PFA (\it{g}\rm)', 'FontSize', 19, 'Color', 'k', 'Interpreter', 'tex');


ylabel('Story loss (USD)', 'FontSize', 19, 'Color', 'k');
% title('SLFs for Column Base Plates', 'FontSize', 14);

% % Enhance the plot box
% set(gca, 'LineWidth', 1.5); % Make the box thicker
% set(gca, 'FontSize', 12); % Setting font size for axes
% set(gca, 'Box', 'on'); % Ensure the box is on
% set(gca, 'LineWidth', 1.2, 'Box', 'on', 'XColor', 'k', 'YColor', 'k');
% Adjust tick labels font size and axes properties
set(gca, 'FontSize', 19, 'LineWidth', 1.5, 'Box', 'on', 'XColor', 'k', 'YColor', 'k');

legend('show', 'Location', 'northwest', 'FontSize', 15);

% Grid on for better readability
% grid on;

% % Set the axes for better visualization
% axis tight;
% ylim([0 1130000]); % Adjust if necessary
% Set the axes for better visualization
axis tight;
ylim([0 1550000]); % Adjust if necessary

% yticks([0, 250000, 500000, 750000, 1000000]);
xlim([0 4]);
xticks([0, 1, 2, 3, 4]);
% Save the figure
% saveas(gcf, 'SLF_vs_IDR_NSA.png'); % Save as PNG
print(gcf, 'SLF_vs_IDR_NSA', '-dpng', '-r1200');
