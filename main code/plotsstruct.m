% Script 3: Band Structure Plotting
% Visualizes the computed band structure
clear; clc;
load('bandstructure.mat'); % Load computed data
load('init_params.mat'); % Load parameters

% Calculate cumulative distance along k-path for x-axis
k_dist = zeros(size(k_points, 1), 1);
for i = 2:size(k_points, 1)
    k_dist(i) = k_dist(i-1) + norm(k_points(i,:) - k_points(i-1,:));
end

% Find positions of high symmetry points
segment_points = [1];
points_per_segment = num_kpoints;
for i = 1:size(BZ_path,1)-1
    segment_points = [segment_points, i*points_per_segment];
end

% Prepare plot
figure('Position', [100, 100, 800, 600]); % Larger figure size
hold on;
box on;

% Plot vertical lines at symmetry points
for i = 1:length(segment_points)
    xline(k_dist(segment_points(i)), '--k', 'LineWidth', 0.5, 'Alpha', 0.5);
end

% Plot bands
for i = 1:size(energy_bands, 2)
    plot(k_dist, energy_bands(:,i), 'LineWidth', 2);
end

% Adjust plot appearance
xlabel('Wave Vector', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy (eV)', 'FontSize', 12, 'FontWeight', 'bold');
title('Electronic Band Structure', 'FontSize', 14, 'FontWeight', 'bold');

% Set symmetry point labels
BZ_labels = {'\Gamma', 'X', 'W', '\Gamma', 'L'};
xticks(k_dist(segment_points));
xticklabels(BZ_labels);

% Add reference line at E = 0
yline(0, '--k', 'LineWidth', 0.5, 'Alpha', 0.5);

% Adjust axis limits for better visualization
ylim([-15, 15]); % Adjust based on your energy range
grid on;

% Create legend
legend(arrayfun(@(x) sprintf('Band %d', x), 1:size(energy_bands, 2), ...
    'UniformOutput', false), 'Location', 'bestoutside');

% Enhance plot appearance
set(gca, 'FontSize', 10);
set(gca, 'LineWidth', 1.5);

% Save the figure
saveas(gcf, 'band_structure.png');
saveas(gcf, 'band_structure.fig');

disp('Plotting complete. Figure saved as band_structure.png and band_structure.fig');