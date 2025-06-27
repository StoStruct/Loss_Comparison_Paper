clc;
clear;
close all;

% Load the data
file1 = 'CombinedOutput_SLF_NSD.csv';
NSD = readmatrix(file1);

% Extract columns for IDR and SLFs
IDR = NSD(:,1);
SLF_1 = NSD(:,2);
SLF_2 = NSD(:,3);
SLF_3 = NSD(:,4);

% Set the font to Times New Roman
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');

% Create figure
figure('Position', [100, 100, 400, 350]);
hold on;

% Plotting each SLF against IDR
plot(IDR, SLF_1, 'LineWidth', 2, 'DisplayName', 'BW-20');
plot(IDR, SLF_2, 'LineWidth', 2, 'DisplayName', 'BW-30', 'LineStyle', '--');
plot(IDR, SLF_3, 'LineWidth', 2, 'DisplayName', 'BW-40', 'LineStyle', ':');

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

% % Set the axes for better visualization
% axis tight;
% ylim([0 1130000]); % Adjust if necessary
% Set the axes for better visualization
axis tight;
ylim([0 1000000]); % Adjust if necessary

yticks([0, 250000, 500000, 750000, 1000000]);
xlim([0 0.2]);
xticks([0, 0.05, 0.1, 0.15, 0.2]);  % Set y-axis ticks HERE 0.2 MEANS RATIO and 20 percent

% Save the figure
% saveas(gcf, 'SLF_vs_IDR_NSD.png'); % Save as PNG
print(gcf, 'SLF_vs_IDR_NSD', '-dpng', '-r1200');
