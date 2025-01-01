% Constants
hbar = 1.0545718e-34; % Reduced Planck constant in J*s
m = 9.10938356e-31;   % Electron mass in kg
a = 1e-10;            % Lattice constant in meters (example value)

% 1. Kinetic Energy in 2D Brillouin Zone
% Coordinates in 2D BZ for a square lattice
kx_corner = pi / a;
ky_corner = pi / a;
kx_side = pi / a;
ky_side = 0;

% Kinetic Energy Calculations
E_corner = (hbar^2 / (2 * m)) * (kx_corner^2 + ky_corner^2);
E_side = (hbar^2 / (2 * m)) * (kx_side^2 + ky_side^2);

% Displaying Results
disp('Kinetic Energy in 2D Brillouin Zone:')
fprintf('E at corner of BZ: %.3e J\n', E_corner);
fprintf('E at side of BZ: %.3e J\n', E_side);
fprintf('Energy ratio (corner/side): %.2f\n', E_corner / E_side);

% 2. Kinetic Energy Ratio for 3D Brillouin Zone
% Coordinates in 3D BZ for a cubic lattice
kx_corner3D = pi / a;
ky_corner3D = pi / a;
kz_corner3D = pi / a;
kx_face3D = pi / a;
ky_face3D = pi / a;
kz_face3D = 0;

% Kinetic Energy Calculations in 3D BZ
E_corner3D = (hbar^2 / (2 * m)) * (kx_corner3D^2 + ky_corner3D^2 + kz_corner3D^2);
E_face3D = (hbar^2 / (2 * m)) * (kx_face3D^2 + ky_face3D^2 + kz_face3D^2);

% Displaying Results
disp('Kinetic Energy in 3D Brillouin Zone:')
fprintf('E at corner of BZ (3D): %.3e J\n', E_corner3D);
fprintf('E at face center of BZ (3D): %.3e J\n', E_face3D);
fprintf('Energy ratio (corner/face center): %.2f\n', E_corner3D / E_face3D);

% 3. Plotting Dispersion Relations for 2D BZ
k2_vals = linspace(-pi/a, pi/a, 100); % k2 values for plotting
k1_zero = 0;
k1_pi = pi / a;

% Energy as a function of k2 with k1 = 0
E_k1_zero = (hbar^2 / (2 * m)) * (k1_zero^2 + k2_vals.^2);

% Energy as a function of k2 with k1 = pi/a
E_k1_pi = (hbar^2 / (2 * m)) * (k1_pi^2 + k2_vals.^2);

% Plotting the results
figure;
plot(k2_vals, E_k1_zero, 'b', 'LineWidth', 1.5); hold on;
plot(k2_vals, E_k1_pi, 'r', 'LineWidth', 1.5);
xlabel('k_2 (1/m)');
ylabel('Energy E (J)');
legend('k_1 = 0', 'k_1 = \pi/a');
title('Dispersion Relations in 2D Brillouin Zone');
grid on;

% Indicate that calculations are complete
disp('Plots generated for dispersion relations in 2D Brillouin Zone.');

% Define lattice basis vectors (in units of lattice spacing 'a')
lattice_type = 'fcc'; % Options: 'simple_cubic', 'bcc', 'fcc', 'hexagonal'
a = 1;  % Lattice spacing
max_order = 2;  % Maximal order of neighboring lattice nodes to be kept

% Basis vectors for different lattice types
switch lattice_type
    case 'simple_cubic'
        basis_vectors = [1 0 0; 0 1 0; 0 0 1] * a;
    case 'bcc'
        basis_vectors = [1 0 0; 0 1 0; 0 0 1; 1/2 1/2 1/2] * a;
    case 'fcc'
        basis_vectors = [1 0 0; 0 1 0; 0 0 1; 1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2] * a;
    case 'hexagonal'
        basis_vectors = [1 0 0; -1/2 sqrt(3)/2 0; 0 0 1.63] * a; % c/a ratio for hexagonal
end

% Generate lattice points by combinations of basis vectors up to max_order
lattice_points = [];
for nx = -max_order:max_order
    for ny = -max_order:max_order
        for nz = -max_order:max_order
            % Calculate lattice vector for each combination of (nx, ny, nz)
            R = nx * basis_vectors(1, :) + ny * basis_vectors(2, :) + nz * basis_vectors(3, :);
            lattice_points = [lattice_points; R, norm(R)]; % Append vector and its magnitude
        end
    end
end

% Sort lattice points by distance from origin
lattice_points = sortrows(lattice_points, 4);

% Display sorted lattice vectors
disp('Sorted lattice vectors by increasing magnitude:')
disp(array2table(lattice_points, 'VariableNames', {'Rx', 'Ry', 'Rz', 'Magnitude'}))

% Plot lattice points in 3D (optional)
figure;
scatter3(lattice_points(:, 1), lattice_points(:, 2), lattice_points(:, 3), 50, lattice_points(:, 4), 'filled');
colorbar;
title('Lattice Points Colored by Distance from Origin');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
axis equal;
