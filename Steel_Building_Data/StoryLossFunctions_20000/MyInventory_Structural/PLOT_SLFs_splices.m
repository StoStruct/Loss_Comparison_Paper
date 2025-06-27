clc;
clear;
close all;

% Load the data
file1 = 'MyInventory_Structural_Splices\CombinedOutput_SLF_Struct_Splices.csv';
ColBasePlate = readmatrix(file1);

% Extract columns for IDR and SLFs
IDR = ColBasePlate(:,1);
SLF_1 = ColBasePlate(:,2);
SLF_2 = ColBasePlate(:,3);
SLF_3 = ColBasePlate(:,4);
SLF_4 = ColBasePlate(:,5);
SLF_5 = ColBasePlate(:,6);
SLF_6 = ColBasePlate(:,7);
SLF_7 = ColBasePlate(:,8);
SLF_8 = ColBasePlate(:,9);
SLF_9 = ColBasePlate(:,10);

% Set the font to Times New Roman
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');

% Create figure
figure('Position', [100, 100, 400, 350]);
hold on;

% Plotting each SLF against IDR
plot(IDR, SLF_1, 'LineWidth', 2, 'DisplayName', 'Splices-36');
% plot(IDR, SLF_2, 'LineWidth', 2, 'DisplayName', 'SLF 2', 'LineStyle', '--');
% plot(IDR, SLF_3, 'LineWidth', 2, 'DisplayName', 'SLF 3', 'LineStyle', ':');
% plot(IDR, SLF_4, 'LineWidth', 2, 'DisplayName', 'SLF 4', 'LineStyle', '-.');
% plot(IDR, SLF_5, 'LineWidth', 2, 'DisplayName', 'SLF 5', 'LineStyle', '-');
% plot(IDR, SLF_6, 'LineWidth', 2, 'DisplayName', 'SLF 6', 'LineStyle', '--');
% plot(IDR, SLF_7, 'LineWidth', 2, 'DisplayName', 'SLF 7', 'LineStyle', ':');
% plot(IDR, SLF_8, 'LineWidth', 2, 'DisplayName', 'SLF 8', 'LineStyle', '-.');
% plot(IDR, SLF_9, 'LineWidth', 2, 'DisplayName', 'SLF 9', 'LineStyle', '-');

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

% legend('show', 'Location', 'northwest', 'FontSize', 15);

% Grid on for better readability
% grid on;

% Set the axes for better visualization
axis tight;
ylim([0 2000000]); % Adjust if necessary
xlim([0 0.2]);
xticks([0, 0.05, 0.1, 0.15, 0.2]);  % Set y-axis ticks
% xticks([0, 0.1, 0.2, 0.3, 0.4]);  % Set y-axis ticks
% Save the figure
print(gcf, 'SLF_vs_IDR_ColSplices', '-dpng', '-r1200');
