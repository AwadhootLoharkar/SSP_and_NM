% Parameters
a = 4.08e-10; % Lattice spacing in meters (4.08 Angstroms for Au)
hbar = 1.0545718e-34; % Reduced Planck's constant (Joule*seconds)
me = 9.10938356e-31; % Electron mass (kg)
eV_conversion = 1.60218e-19; % Joules to eV conversion factor

% Define a larger range of k values (e.g., from -6 to 6 in units of 2π/a)
k_vals = linspace(-6 * pi / a, 6 * pi / a, 2000); % More points for smoothness

% Define number of bands to plot (for each reciprocal lattice point)
num_bands = 4; % Number of bands to display
E_bands = zeros(num_bands, length(k_vals));

% Calculate the energy for each band (λ = 1, 2, 3, 4) over the extended range
for n = 1:num_bands
    for j = 1:length(k_vals)
        k = k_vals(j);
        % Calculate energy for each band as a parabolic function offset by (2πn/a)
        k_shift = (n - 1) * (2 * pi / a); % Shifted wave vector for each band
        E_bands(n, j) = hbar^2 * (k - k_shift)^2 / (2 * me) / eV_conversion; % Convert to eV
    end
end

% Plotting
figure;
hold on;

% Plot each band with a different line style and color
plot(k_vals / (2 * pi / a), E_bands(1, :), 'r-', 'LineWidth', 1.5); % Band 1 (solid red)
plot(k_vals / (2 * pi / a), E_bands(2, :), 'b--', 'LineWidth', 1.5); % Band 2 (dashed blue)
plot(k_vals / (2 * pi / a), E_bands(3, :), 'g:', 'LineWidth', 1.5); % Band 3 (dotted green)
plot(k_vals / (2 * pi / a), E_bands(4, :), 'c-.', 'LineWidth', 1.5); % Band 4 (dash-dot light blue)

% Add plot details
xlabel('k (Units 2\pi/a)');
ylabel('E(k) (eV)');
title('Extended Energy Bands in a One-Dimensional Empty Lattice Model');
ylim([0, 35]); % Adjust y-axis range to match the example image
xlim([-6, 6]); % Larger k range to show multiple BZs
grid on;

% Display
hold off;
