% Step 1: Define appropriate units
recipunit = 1e10;   % 1 / m to Å^-1 conversion factor
hbar = 1.055e-34;   % Reduced Planck's constant in J.s
m_e = 9.11e-31;     % Mass of an electron in kg
ekinscale = hbar^2 / (2 * m_e); % Scale factor to convert k^2 to energy in Joules
E0 = 1.6e-19;      % Energy of electron (J)

% Calculate the wave vector k
k = sqrt((2 * m_e * E0) / (hbar^2));

% Convert ekinscale to eV (1 eV = 1.60218e-19 J)
ekinscale = ekinscale / 1.60218e-19; % Scale factor in eV

% Step 2: Define grid of points and energy bounds

% Define the limits of the barrier
x_min = 0;  % Minimum x' (Å)
x_max = 5;  % Maximum x' (Å)
dx = 0.1;   % Grid spacing (Å)

% Create a grid of points inside the barrier
x_grid = x_min:dx:x_max;  % Grid points

% Define energy bounds
E_a = 0;    % Lower bound of energy (eV)
E_b = 3;    % Upper bound of energy (eV)

% Display the grid and energy bounds
disp('Grid points:');
disp(x_grid);
disp(['Energy bounds: Ea = ', num2str(E_a), ' eV, Eb = ', num2str(E_b), ' eV']);



% Step 3: Define the spatial part of the incident plane wave

% Define the amplitude (arbitrarily chosen, can be 1 for simplicity)
A = 1;  

% Calculate the wave function for the incident plane wave
Phi_0_k = A * exp(1i * k * x_grid);  % Using complex exponential

% Display results
fprintf('Wave vector k: %.4e 1/Angstrom\n', k);
fprintf('Incident plane wave function (first 5 points):\n');
disp(Phi_0_k(1:5));

% Step 4: Define the perturbation by setting U0 and computing V(x)

% Define the height of the potential barrier (in eV)
U0 = 2;  % Height of the barrier

% Initialize the potential array V(x)
V = zeros(size(x_grid));  % Start with V(x) = 0 everywhere

% Set the potential inside the barrier
V(x_grid >= x_min & x_grid <= x_max) = U0; 

% Display results
fprintf('Potential V(x) (first 10 points):\n');
disp(V(1:10));  % Display the first 10 points of V(x)

% Step 5: Define the Green's function for all pairs of points inside the perturbation

% Initialize the Green's function matrix
G = zeros(length(x_grid), length(x_grid));

% Define the wave number for free space (outside the barrier)
k_free = sqrt(2 * m_e * E0 / hbar^2);  % Wave number for E0

% Define the wave number inside the barrier (complex if E < U0)
k_inside = sqrt(2 * m_e * (E0 - U0) / hbar^2);  % Wave number for E < U0

% Loop through all combinations of x and x' to compute G(x, x')
for i = 1:length(x_grid)
    for j = 1:length(x_grid)
        if x_grid(i) < x_min || x_grid(i) > x_max  % Outside the barrier
            % Green's function in free space
            G(i,j) = (1/(2*k_free)) * exp(1i*k_free*abs(x_grid(i)-x_grid(j)));
        elseif x_grid(i) >= x_min && x_grid(i) <= x_max  % Inside the barrier
            % Green's function in the barrier
            G(i,j) = (1/(2*k_inside)) * exp(1i*k_inside*abs(x_grid(i)-x_grid(j)));
        end
    end
end

% Display results
fprintf('Greens function G(x, x_prime) (first 10x10 block):\n');
disp(G(1:10, 1:10));  % Display the first 10x10 block of G(x, x')

% Step 6: Build the scattering matrix

% Define the number of points in the grid
n_points = length(x_grid); % Number of grid points

% Initialize the scattering matrix S with zeros
S = zeros(n_points); 

% Loop through all combinations of grid points to populate the scattering matrix
for i = 1:n_points
    for j = 1:n_points
        % Compute the scattering matrix using the Green's function
        % This is a simplistic representation; the actual relation may vary
        S(i, j) = G(i, j) * exp(-1i * k * abs(x_grid(i) - x_grid(j))); 
    end
end

% Display the scattering matrix
fprintf('Scattering matrix S (first 10x10 block):\n');
disp(S(1:10, 1:10));  % Display the first 10x10 block of the scattering matrix

% Step 7: Solve for the scattering eigenstate wavefunction Φ_k(x_i)

% Initialize the eigenstate wavefunction array
Phi_k = zeros(size(x_grid)); 

% Calculate the eigenstate wavefunction using the scattering matrix
for i = 1:n_points
    Phi_k(i) = sum(S(i, :) .* Phi_0_k); % Sum over all states with incident wavefunction
end

% Display results
fprintf('Scattering eigenstate wavefunction (first 10 points):\n');
disp(Phi_k(1:10));  % Display the first 10 points of the eigenstate wavefunction

% Step 8: Define a larger grid of points and compute the wavefunction

% Define a larger grid for x spanning from -x_max to x_max (including space outside the barrier)
x_large_grid = linspace(-x_max, x_max, 100);  % 100 points from -x_max to x_max

% Initialize the wavefunction array for the larger grid
Phi_k_large = zeros(size(x_large_grid));

% Compute the wavefunction at all points in the larger grid
for i = 1:length(x_large_grid)
    if x_large_grid(i) >= x_min && x_large_grid(i) <= x_max
        % Find the closest point in the original grid to use the scattering matrix S
        [~, idx] = min(abs(x_large_grid(i) - x_grid));  % Get the index of the closest point in x_grid
        Phi_k_large(i) = sum(S(idx, :) .* Phi_0_k);  % Use the closest point in the original grid
    else
        % Outside the barrier (free space)
        Phi_k_large(i) = exp(1i * k * x_large_grid(i));  % Free space wavefunction
    end
end

% Display the wavefunction on the larger grid (first 10 points)
fprintf('Wavefunction on the larger grid (first 10 points):\n');
disp(Phi_k_large(1:10));

% Step 9: Visualize the wavefunction

% Plot the wavefunction on the larger grid
figure;
plot(x_large_grid, real(Phi_k_large), 'r', 'LineWidth', 1.5);  % Plot the real part
hold on;
plot(x_large_grid, imag(Phi_k_large), 'b', 'LineWidth', 1.5);  % Plot the imaginary part
title('Wavefunction \Phi_k(x) on the Larger Grid');
xlabel('x (Å)');
ylabel('\Phi_k(x)');
legend('Real part', 'Imaginary part');
grid on;

% Step 10: Calculate the transmission and reflection coefficients

% Extract the transmitted part of the wavefunction (far right of the barrier)
Phi_transmitted = Phi_k_large(x_large_grid > x_max);

% Extract the reflected part of the wavefunction (far left of the barrier)
Phi_reflected = Phi_k_large(x_large_grid < x_min);

% Calculate the transmission coefficient (squared modulus of transmitted wave)
T = abs(Phi_transmitted(end))^2;

% Calculate the reflection coefficient (squared modulus of reflected wave)
R = abs(Phi_reflected(1))^2;

% Display results
fprintf('Transmission coefficient (T): %.4f\n', T);
fprintf('Reflection coefficient (R): %.4f\n', R);

% Step 11: Compute the Reflection and Transmission Coefficients

% Reflection coefficient (R): Compare the wavefunction in the negative x region (-x_max)
% with the incident plane wave (since it is reflected back)
reflected_wave = Phi_k_large(x_large_grid < 0);
incident_wave = exp(1i * k * x_large_grid(x_large_grid < 0));  % Free space incident wave

R = abs(sum(reflected_wave) / sum(incident_wave))^2;  % Reflection coefficient (magnitude squared)

% Transmission coefficient (T): Compare the wavefunction in the positive x region (+x_max)
% with the incident plane wave (since it is transmitted through the barrier)
transmitted_wave = Phi_k_large(x_large_grid > 0);
T = abs(sum(transmitted_wave) / sum(incident_wave))^2;  % Transmission coefficient (magnitude squared)

% Display reflection and transmission coefficients
fprintf('Reflection coefficient R: %.4f\n', R);
fprintf('Transmission coefficient T: %.4f\n', T);

% Step 12: Normalize the Wavefunction

% Compute the norm (integral of |Phi_k_large|^2 over all space)
norm_factor = sqrt(sum(abs(Phi_k_large).^2) * (x_large_grid(2) - x_large_grid(1)));  % Sum with spacing dx

% Normalize the wavefunction
Phi_k_large_normalized = Phi_k_large / norm_factor;

% Display normalized wavefunction (first 10 points)
fprintf('Normalized wavefunction (first 10 points):\n');
disp(Phi_k_large_normalized(1:10));

% Step 13: Plot the Potential and the Wavefunction

% Plot the potential V(x) across the entire grid
figure;
plot(x_large_grid, [zeros(1, length(x_large_grid(x_large_grid < x_min))) V zeros(1, length(x_large_grid(x_large_grid > x_max)))], 'r', 'LineWidth', 2);
hold on;

% Plot the real part of the wavefunction
plot(x_large_grid, real(Phi_k_large_normalized), 'b', 'LineWidth', 2);
hold on;

% Plot the imaginary part of the wavefunction
plot(x_large_grid, imag(Phi_k_large_normalized), 'g', 'LineWidth', 2);

% Labels and title
xlabel('x (Å)');
ylabel('Amplitude');
title('Potential Barrier and Normalized Wavefunction');
legend('Potential V(x)', 'Re(\Phi_k)', 'Im(\Phi_k)');
grid on;

% Display the plot
hold off;
