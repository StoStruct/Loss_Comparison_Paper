clc;
clear;
close all;

% Load the data
file1 = 'MyInventory_Structural\\MyInventory_Structural_Splices\CombinedOutput_SLF_Struct_ColBasePlate.csv';
ColBasePlate = readmatrix(file1);

% Extract columns for IDR and SLFs
IDR = ColBasePlate(:,1);
SLF_1 = ColBasePlate(:,2);
SLF_2 = ColBasePlate(:,5);
SLF_3 = ColBasePlate(:,8);

% Create a high-quality figure
figure;
hold on; % Hold on to plot multiple lines

% Plotting each SLF against IDR
plot(IDR, SLF_1, 'LineWidth', 2, 'DisplayName', 'SLF 1');
plot(IDR, SLF_2, 'LineWidth', 2, 'DisplayName', 'SLF 2', 'LineStyle', '--');
plot(IDR, SLF_3, 'LineWidth', 2, 'DisplayName', 'SLF 3', 'LineStyle', ':');

% Adding labels and title
xlabel('Inter-Storey Drift Ratio (IDR)', 'FontSize', 12);
ylabel('Storey Loss Function (SLF)', 'FontSize', 12);
title('SLFs for Column Base Plates', 'FontSize', 14);

% Enhance the plot box
set(gca, 'LineWidth', 1.2); % Make the box thicker
set(gca, 'FontSize', 12); % Setting font size for axes
set(gca, 'Box', 'on'); % Ensure the box is on

% Adjust legend
legend('show', 'Location', 'best');

% Grid on for better readability
% grid on;

% Set the axes for better visualization
axis tight;
ylim([0 1000000]); % Adjust if necessary

% Save the figure
saveas(gcf, 'MyInventory_Structural\SLF_vs_IDR_ColBasePlate.png'); % Save as PNG
