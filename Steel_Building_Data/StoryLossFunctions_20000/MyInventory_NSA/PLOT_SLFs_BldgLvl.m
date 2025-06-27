clc;
clear;
close all;

% Load the data
file1 = 'MyInventory_NSA_BuildingLevelComponents\CombinedOutput_SLF_NSA_Bldglvl.csv';
BldgLvl = readmatrix(file1);

% Extract columns for IDR and SLFs
IDR = BldgLvl(:,1);
SLF_1 = BldgLvl(:,2);
SLF_2 = BldgLvl(:,3);
SLF_3 = BldgLvl(:,4);
SLF_4 = BldgLvl(:,5);
SLF_5 = BldgLvl(:,6);
SLF_6 = BldgLvl(:,7);
SLF_7 = BldgLvl(:,8);
SLF_8 = BldgLvl(:,9);
SLF_9 = BldgLvl(:,10);
SLF_10 = BldgLvl(:,11);
SLF_11 = BldgLvl(:,12);
SLF_12 = BldgLvl(:,13);
SLF_13 = BldgLvl(:,14);
SLF_14 = BldgLvl(:,15);
SLF_15 = BldgLvl(:,16);

% Set the font to Times New Roman
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');

% Create figure
% figure('Position', [100, 100, 400, 350]);
figure('Position', [100, 100, 390, 350]);
hold on;

% % Plotting each SLF against IDR
% plot(IDR, SLF_1, 'LineWidth', 2, 'DisplayName', 'SC-1;BW-20');
% % plot(IDR, SLF_2, 'LineWidth', 2, 'DisplayName', 'SC-1, BW-30', 'LineStyle', '--');
% % plot(IDR, SLF_3, 'LineWidth', 2, 'DisplayName', 'SC-1, BW-40', 'LineStyle', ':');
% plot(IDR, SLF_4, 'LineWidth', 2, 'DisplayName', 'SC-5, BW-20', 'LineStyle', '-.');
% % plot(IDR, SLF_5, 'LineWidth', 2, 'DisplayName', 'SC-5, BW-30', 'LineStyle', '-');
% % plot(IDR, SLF_6, 'LineWidth', 2, 'DisplayName', 'SC-5, BW-40', 'LineStyle', '--');
% plot(IDR, SLF_7, 'LineWidth', 2, 'DisplayName', 'SC-9, BW-20', 'LineStyle', ':');
% % plot(IDR, SLF_8, 'LineWidth', 2, 'DisplayName', 'SC-9, BW-30', 'LineStyle', '-.');
% % plot(IDR, SLF_9, 'LineWidth', 2, 'DisplayName', 'SC-9, BW-40', 'LineStyle', '-');
% plot(IDR, SLF_10, 'LineWidth', 2, 'DisplayName', 'SC-14, BW-20', 'LineStyle', ':');
% % plot(IDR, SLF_11, 'LineWidth', 2, 'DisplayName', 'SC-14, BW-30', 'LineStyle', '-.');
% % plot(IDR, SLF_12, 'LineWidth', 2, 'DisplayName', 'SC-14, BW-40', 'LineStyle', '-');
% plot(IDR, SLF_13, 'LineWidth', 2, 'DisplayName', 'SC-19, BW-20', 'LineStyle', ':');
% % plot(IDR, SLF_14, 'LineWidth', 2, 'DisplayName', 'SC-19, BW-30', 'LineStyle', '-.');
% % plot(IDR, SLF_15, 'LineWidth', 2, 'DisplayName', 'SC-19, BW-40', 'LineStyle', '-');

% Define distinct colors (ColorBrewer-style palette)
colorSet = [
    0.1216, 0.4667, 0.7059;   % blue
    0.8902, 0.4667, 0.7608;   % pink
    0.1725, 0.6275, 0.1725;   % green
    1.0000, 0.4980, 0.0549;   % orange
    0.8392, 0.1529, 0.1569;   % red
];

% Distinct linestyle patterns
styleSet = {'-', '-.', ':', '--', '-'};  % Customize more if needed

hold on
plot(IDR, SLF_2,  'LineWidth', 2, 'Color', colorSet(1,:), 'LineStyle', styleSet{1}, 'DisplayName', 'SC-1; BW-30');
plot(IDR, SLF_5,  'LineWidth', 2, 'Color', colorSet(2,:), 'LineStyle', styleSet{2}, 'DisplayName', 'SC-5; BW-30');
plot(IDR, SLF_8,  'LineWidth', 2, 'Color', colorSet(3,:), 'LineStyle', styleSet{3}, 'DisplayName', 'SC-9; BW-30');
plot(IDR, SLF_11, 'LineWidth', 2, 'Color', colorSet(4,:), 'LineStyle', styleSet{4}, 'DisplayName', 'SC-14; BW-30');
plot(IDR, SLF_14, 'LineWidth', 2, 'Color', colorSet(5,:), 'LineStyle', styleSet{5}, 'DisplayName', 'SC-19; BW-30');

% legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10);


% Adding labels and title
% xlabel('PFA', 'FontSize', 19, 'Color', 'k');
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

% legend('show', 'Location', 'northoutside', 'FontSize', 15);

% Grid on for better readability
% grid on;

% % Set the axes for better visualization
% axis tight;
% ylim([0 1130000]); % Adjust if necessary
% Set the axes for better visualization
axis tight;
ylim([0 1500000]); % Adjust if necessary
% yticks([0, 250000, 500000, 750000, 1000000]);
% yticks([0, 250000, 500000, 750000, 1000000]);
xlim([0 4]);
xticks([0, 1, 2, 3, 4]);  % Set y-axis ticks

% Save the figure
% saveas(gcf, 'SLF_vs_IDR_BldgLvl.png'); % Save as PNG
print(gcf, 'SLF_vs_IDR_BldgLvl', '-dpng', '-r1200');
